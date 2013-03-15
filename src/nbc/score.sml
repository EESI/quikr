signature SCORE = sig
	val score: int * real * (string -> int option) * real * string -> real
end

structure Score :> SCORE = struct
	fun |> (x, f) = f x
	infix |>

	fun addCount (hitsum, fcount, gcount, totalWords) =
		Kahan.add (hitsum, Real.fromInt fcount * Math.ln (Real.fromInt gcount / totalWords))
	fun addNmer (totalWords, getGenomeCount) (nmer, ref fcount, (misses, anyhits, hitsum)) =
		case getGenomeCount nmer of
			NONE => (misses + 1, anyhits, hitsum)
			| SOME gcount => (
				misses, true
				, addCount (hitsum, fcount, gcount, totalWords)
			)
	fun score (order, missConstant, getGenomeCount, totalWords, fragment) =
		let
			val add = addNmer (totalWords, getGenomeCount)
			val seed = (0, false, Kahan.zero)
			val (misses, anyhits, hitsum) =
				Nmer.count (order, fragment) |> HashTable.foldi add seed
		in
			if anyhits then
				Kahan.add (hitsum, Math.ln missConstant * Real.fromInt misses)
				|> Kahan.sum
			else Real.negInf
		end
end
