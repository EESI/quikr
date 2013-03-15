structure Main = struct

fun |> (x, f) = f x
infix |>

val order = Options.order
val genomes = case Options.genomes of
	SOME x => x |> TextIO.openIn |> Misc.sequenceLines
	| NONE => Options.genomesDir |> Misc.sortedDirectoryNoPrefix |> Sequence.fromList
fun input () =
	case Fasta.sequence (Gzip.openIn Options.input) of
		NONE => Fail.fail "input file is not FASTA format"
		| SOME x =>
			Sequence.map (fn (header, data) =>
				(header, String.map Char.toUpper data)
			) x
fun output genome =
	let
		fun inner format = case format of
			Options.Matlab {variable, file} =>
				let
					val matlab = Matlab.openOut file
					val doubleArray = Matlab.beginDoubleArray (matlab, variable)
				in {
					write = fn (_, score) =>
						Matlab.writeDouble (doubleArray, score)
					, close = fn () => (
						Matlab.concludeDoubleArray doubleArray
						; Matlab.closeOut matlab
					)
				} end
			| Options.Text filename =>
				let
					val gzip = Gzip.openOut filename
				in {
					write = fn (header, score) =>
						TextIO.output (
							gzip
							, Real.fmt (StringCvt.FIX (SOME 8)) score
						)
					, close = fn () => TextIO.closeOut gzip
				} end
			| Options.Dual {text, matlab} =>
				let
					val text = inner (Options.Text text)
					val matlab = inner (Options.Matlab matlab)
				in {
					write = fn x => (
						#write text x
						; #write matlab x
					), close = fn () => (
						#close text ()
						; #close matlab ()
					)
				} end
	in
		inner (Options.output genome)
	end
fun totalWords (genome, order) =
        let
		val name = Options.totalWords (genome, order)
		val input = TextIO.openIn name
	in
		(case Option.mapPartial Real.fromString (Misc.inputLine input) of
			NONE => Fail.fail ("could not read number from " ^ name)
			| SOME r => r
		) before TextIO.closeIn input
	end

val () =
	Sequence.app (fn gname =>
		let
			val {write, close} = output gname
			val totalWords = totalWords (gname, order)
			val genome = (
				Stopwatch.start ("Loading genome " ^ gname)
				; Genome.load (gname, order)
				before Stopwatch.finish ()
			)
		in
			Stopwatch.start ("Scoring fragments")
			; input () |> Sequence.app (fn (header, fragment) =>
				write (
					header
					, Score.score (
						order
						, Options.missConstant
						, fn nmer => Genome.get (genome, nmer)
						, totalWords
						, fragment
					)
				)
			); Stopwatch.finish ()
			; close ()
		end
	) genomes

end
