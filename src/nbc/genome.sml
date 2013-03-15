signature GENOME = sig
	exception Bad
	type t
	val load: string * int -> t
	val get: t * string -> int option
end

structure Genome :> GENOME = struct
	exception Bad
	fun |> (x, f) = f x
	infix |>

	type t = (string, int) HashTable.hash_table
	fun load (gname, order) =
		let
			val h = HashTable.mkTable
				(HashString.hashString, op =)
				(1024 * 1024, Fail "")
		in
			Options.genomeText (order, gname) |> Gzip.openIn |> Misc.sequenceLines
			|> Sequence.map (fn s => (
				 case Misc.split2 s of
					SOME (count, nmer) => (
						nmer
						, case Int.fromString count of
							NONE => raise Bad
							| SOME x => x
					) | NONE => raise Bad
			)) |> Sequence.app (HashTable.insert h)
			; h
		end
	fun get (h, nmer) = HashTable.find h nmer
end
