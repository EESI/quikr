structure History :> sig
	type history
	val create: int -> history
	val clear: history -> history
	val push: history * char -> string option * history
end = struct
	type history = {maximumSize: int, content: string}
	fun create maximumSize = {maximumSize = maximumSize, content = ""}
	fun clear {maximumSize, content = _} =
		{maximumSize = maximumSize, content = ""}
	fun addToString (string, newCharacter) =
		let
			val size = size string
		in
			CharVector.tabulate (
				size + 1
				, fn index =>
					if index = size then newCharacter
					else String.sub (string, index)
			)
		end
	fun shiftIntoString (string, newCharacter) =
		let
			val size = size string
		in
			CharVector.tabulate (
				size
				, fn index =>
					if index = size - 1 then newCharacter
					else String.sub (string, index + 1)
			)
		end
	fun push ({maximumSize, content}, newCharacter) =
		let
			val currentSize = size content
		in
			if currentSize = maximumSize then
				let
					val newContent = shiftIntoString (
						content
						, newCharacter
					)
				in (
					SOME newContent
					, {
						maximumSize = maximumSize
						, content = newContent
					}
				) end
			else if currentSize = maximumSize - 1 then
				let
					val newContent = addToString (
						content
						, newCharacter
					)
				in (
					SOME newContent
					, {
						maximumSize = maximumSize
						, content = newContent
					}
				) end
			else (
				NONE
				, {
					maximumSize = maximumSize
					, content = addToString (
						content
						, newCharacter
					)
				}
			)
		end
end
