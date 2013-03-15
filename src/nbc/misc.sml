signature MISC = sig
	val inputLine: TextIO.instream -> string option
	val sequenceFromRead: (TextIO.instream -> 'a option) -> TextIO.instream -> 'a Sequence.t
	val sequenceLines: TextIO.instream -> string Sequence.t
	val sortedDirectoryNoPrefix: string -> string list
	val sortedDirectory: string -> string list
	val substitute: (string * string) list -> string -> string option
	val basenameWithoutExtension: string -> string
	val split2: string -> (string * string) option
	val longestCommonSubstring: string list -> string
end

structure Misc :> MISC = struct
	fun |> (x, f) = f x
	infix |>
	fun \ f x y = f (x, y)
	fun sequenceFromRead f ioc = Sequence.from (fn () => f ioc)
	fun inputLine instream = case TextIO.inputLine instream of
		NONE => NONE
		| SOME x => SOME (String.substring (x, 0, size x - 1))
	val sequenceLines = sequenceFromRead inputLine
	fun sortedDirectoryNoPrefix d =
		let
			val h = OS.FileSys.openDir d
			fun loop l =
				case OS.FileSys.readDir h of
					NONE => (
						OS.FileSys.closeDir h
						; ListMergeSort.sort (op >) l
					) | SOME n => loop (
						if String.sub (n, 0) = #"." then l
						else n :: l
					)
		in
			loop nil
		end
	fun sortedDirectory d =
		map (\OS.Path.concat d) (sortedDirectoryNoPrefix d)
	fun assoc list key =
		case
			List.find (fn (possibility, _) => possibility = key) list
		of
			NONE => NONE
			| SOME (_, answer) => SOME answer
	local
		exception NotFound
	in
		fun substitute v s =
			Substitution.substitute
				(fn s => (case assoc v s of
					NONE => raise NotFound
					| SOME x => x
				)) s
			handle NotFound => NONE
	end
	fun basenameWithoutExtension s = OS.Path.base (OS.Path.file s)
	local
		fun index f (string, offset) =
			let
				val last = size string
				fun loop i =
					if i = last then ~1
					else if f (String.sub (string, i)) then i
					else loop (i + 1)
			in
				loop offset
			end
	in
		val indexWhitespace = index Char.isSpace
		val indexNonWhitespace = index (not o Char.isSpace)
	end
	fun split2 s =
		let
			val f1b = indexNonWhitespace (s, 0)
		in
			if f1b = ~1 then NONE
			else let
				val f1e = indexWhitespace (s, f1b + 1)
			in
				if f1e = ~1 then NONE
				else let
					val f2b = indexNonWhitespace (s, f1e + 1)
				in
					if f2b = ~1 then NONE
					else let
						val f2e = indexWhitespace (s, f2b + 1)
					in
						if f2e = ~1 then SOME (
							substring (s, f1b, f1e - f1b)
							, substring (s, f2b, size s - f2b)
						) else if indexNonWhitespace (s, f2e + 1) = ~1 then
							SOME (
								substring (s, f1b, f1e - f1b)
								, substring (s, f2b, f2e - f2b)
							)
						else NONE
					end
				end
			end
		end
	fun longestCommonSubstring strings =
		case
			foldl (fn (a, b) => (case b of
				NONE => SOME a
				| SOME b => if size a < size b then SOME a else SOME b
			)) NONE strings
		of
			NONE => ""
			| SOME shortest => let
				val minimumSize = size shortest
				fun loop (size, offset) =
					let
						fun next () =
							if size + offset = minimumSize then
								loop (size - 1, 0)
							else loop (size, offset + 1)
						val sub = substring (shortest, offset, size)
					in
						if List.all (String.isSubstring sub) strings then sub
						else next ()
					end
			in
				loop (minimumSize, 0)
			end
end
