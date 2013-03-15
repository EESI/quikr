signature PARSE_STATE = sig
	type ('argument, 'result) state
	val create: unit -> ('argument, 'result) state
	type ('argument, 'result) handler =
		TextIO.instream * 'argument -> 'result
	datatype whichCharacters = These of char list | Any
	val build: {
		state: ('argument, 'result) state
		, characters:
			(whichCharacters * ('argument, 'result) handler) list
		, endOfFile: ('argument, 'result) handler
	} -> unit
	val enter:
		('argument, 'result) state
			* TextIO.instream
			* 'argument
		-> 'result
end

structure ParseState :> PARSE_STATE = struct
	type ('argument, 'result) handler =
		TextIO.instream * 'argument -> 'result
	datatype whichCharacters = These of char list | Any
	type ('argument, 'result) state = {
		byCharacter: Int8.int vector ref
		, byIndex: ('argument, 'result) handler vector ref
		, endOfFile: ('argument, 'result) handler option ref
	}
	fun create () = {
		byCharacter = ref (Vector.fromList nil)
		, byIndex = ref (Vector.fromList nil)
		, endOfFile = ref NONE
	}
	fun build {
		state = {byCharacter, byIndex, endOfFile}
		, characters
		, endOfFile = newEndOfFile
	} =
		let
			val characters = vector characters
			fun equal (one: char) (two: char) =
				one = two
			fun shallHandle ((whichToHandle, _), char) =
				case whichToHandle of
					Any => true
					| These charactersToHandle =>
						List.exists (equal char)
							charactersToHandle
			fun charToIndex char =
				case
					Vector.findi (fn (_, handler) =>
						shallHandle (handler, char)
					) characters
				of
					NONE => raise Fail (
						"ParseState.build: "
						^ Char.toString char
						^ " not found"
					) | SOME (index, _) =>
						Int8.fromInt index
			fun handlerToFunction (_, function) = function
			fun indexToFunction index = handlerToFunction (
				Vector.sub (characters, index)
			)
		in
			byCharacter := Vector.tabulate (
				Char.maxOrd + 1
				, charToIndex o chr
			); byIndex :=
				Vector.map (fn (_, function) =>
					function
				) characters
			; endOfFile := SOME newEndOfFile
		end
	fun enter (
		{
			byCharacter = ref byCharacter
			, byIndex = ref byIndex
			, endOfFile = ref endOfFile
		}
		, instream
		, argument
	) = case TextIO.input1 instream of
		NONE => (valOf endOfFile) (instream, argument)
		| SOME char => Vector.sub (
			byIndex
			, Int8.toInt (Vector.sub (byCharacter, ord char))
		) (instream, argument)
end
