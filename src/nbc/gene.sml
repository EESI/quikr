signature GENE = sig
	val reverse: string -> string
	val first: int -> string
	val next: string -> string option
end

structure Gene :> GENE = struct
	fun reverse s =
		let
			val n = size s
			val m = n - 1
			fun opposite c = case c of
				#"A" => #"T"
				| #"T" => #"A"
				| #"C" => #"G"
				| #"G" => #"C"
				| _ => c
		in
			CharVector.tabulate (n, fn i => opposite (String.sub (s, m - i)))
		end
	fun first order = CharVector.tabulate (order, fn _ => #"A")
	fun next nmer =
		let
			val order = size nmer
			fun finish (rightmostNonT, replacement) = CharVector.tabulate (
				order
				, fn index =>
					case
						Int.compare (
							index
							, rightmostNonT
						)
					of
						LESS => String.sub (nmer, index)
						| EQUAL => replacement
						| GREATER => #"A"
			)
			fun continue index =
				if index < 0 then NONE
				else case String.sub (nmer, index) of
					#"A" => SOME (finish (index, #"C"))
					| #"C" => SOME (finish (index, #"G"))
					| #"G" => SOME (finish (index, #"T"))
					| #"T" => continue (index - 1)
					| _ => raise Fail "Invalid base"
		in
			continue (size nmer - 1)
		end
end
