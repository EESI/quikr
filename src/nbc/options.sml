signature OPTIONS = sig
	val order: int
	val input: string
	val genomesDir: string
	val genomes: string option
	val genomeText: int * string -> string
	val totalWords: string * int -> string
	val missConstant: real
	datatype output =
		Matlab of {variable: string, file: string}
		| Text of string
		| Dual of {text: string, matlab: {variable: string, file: string}}
	val output: string -> output
	val arguments: string list
end

structure Options :> OPTIONS = struct
	fun |> (x, f) = f x
	infix |>
	datatype outformat = Otext | Omatlab | Odual
	local
		val order = ref 15
		val fasta = ref NONE
		val genomesDir = ref "/var/lib/genomes"
		val genomes = ref NONE
		val genomeText = ref "$genomes_dir/$genome/${order}perword.gz"
		val totalWords = ref "$genomes_dir/$genome/${order}total"
		val missConstant = ref (1.0 / 26067530.0 / 1000000.0)
		val outFormat = ref Otext
		val variable = ref "ps_g"
		val outFile = ref "$input-$order-$genome.$extension"
		fun withDefault (description, default) = concat (
			description
			:: " ("
			:: (case default of
				NONE => ["no default)"]
				| SOME x => ["default: ", x, ")"])
		)
		val options = [
			{
				short = "r", long = ["order"]
				, help = withDefault
					("Set the length of words used in matches", NONE)
				, desc = GetOpt.ReqArg (
					fn x => case Int.fromString x of
						NONE => Fail.fail
							"given order is not a valid integer"
						| SOME y => order := y
					, "integer"
				)
			}, {
				short = "a", long = ["fasta-input"]
				, help = "FASTA input file"
				, desc = GetOpt.ReqArg (fn x => fasta := SOME x, "filename")
			}, {
				short = "j", long = ["genomes-dir"]
				, help = withDefault (
					"Set the directory where genomes are stored"
					, SOME (!genomesDir)
				), desc = GetOpt.ReqArg (fn x => genomesDir := x, "directory")
			}, {
				short = "g", long = ["genome-list"]
				, help = withDefault (
					"Read genomes from this file, one per line"
					, !genomes
				), desc = GetOpt.ReqArg (fn x => genomes := SOME x, "filename")
			}, {
				short = "w", long = ["genome-text"]
				, help = withDefault (
					"Filename to read per-word counts from"
					, SOME (!genomeText)
				), desc = GetOpt.ReqArg (fn x => genomeText := x, "filename")
			}, {
				short = "t", long = ["total-words"]
				, help = withDefault (
					"Filename to read total word counts from"
					, SOME (!totalWords)
				), desc = GetOpt.ReqArg (fn x => totalWords := x, "filename")
			}, {
				short = "k", long = ["miss-constant"]
				, help = withDefault (
					"Set the constant used for scoring missing words"
					, SOME (Real.toString (!missConstant))
				), desc = GetOpt.ReqArg (
					fn x => (case Real.fromString x of 
						SOME y => missConstant := y
						| NONE => Fail.fail
							"given miss constant is not a valid number"
					), "number"
				)
			}, {
				short = "x", long = ["text"]
				, help = "Select gzipped text output (default)"
				, desc = GetOpt.NoArg (fn () => outFormat := Otext)
			}, {
				short = "m", long = ["matlab"]
				, help = "Select Matlab binary output"
				, desc = GetOpt.NoArg (fn () => outFormat := Omatlab)
			}, {
				short = "u", long = ["dual"]
				, help = "Select both gzipped text and Matlab binary output"
				, desc = GetOpt.NoArg (fn () => outFormat := Odual)
			}, {
				short = "v", long = ["variable"]
				, help = withDefault (
					"Set Matlab output variable"
					, SOME (!variable)
				), desc = GetOpt.ReqArg (fn x => variable := x, "name")
			}, {
				short = "o", long = ["output-file"]
				, help = withDefault (
					"Set output filename(s)"
					, SOME (!outFile)
				), desc = GetOpt.ReqArg (fn x => outFile := x, "filename")
			}
		]
		val optionsWithHelp = {
			short = "h", long = ["help"]
			, help = "Display help"
			, desc = GetOpt.NoArg (fn () => (
				print (
					GetOpt.usageInfo {
						header = "score <options>"
						, options = options
					} ^ "\n"
				); OS.Process.exit OS.Process.success
			))
		} :: options
	in
		val (_, arguments) = GetOpt.getOpt {
			argOrder = GetOpt.Permute
			, options = optionsWithHelp
			, errFn = print
		} (CommandLine.arguments ())
		val order = !order
		val input = case !fasta of 
			NONE => Fail.fail "No input file given"
			| SOME x => x
		val genomesDir = !genomesDir
		val genomes = !genomes
		val genomeText = fn (order, genome) => (
			case
				Misc.substitute [
					("genomes_dir", genomesDir)
					, ("genome", genome)
					, ("order", Int.toString order)
				] (!genomeText)
			of
				NONE => Fail.fail "bad per-word count filename syntax"
				| SOME x => x
		)
		val totalWords = fn (genome, order) => (
			case
				Misc.substitute [
					("genomes_dir", genomesDir)
					, ("genome", genome)
					, ("order", Int.toString order)
				] (!totalWords)
			of
				NONE => Fail.fail "bad total word count filename syntax"
				| SOME x => x
		)
		val missConstant = !missConstant
		datatype output =
			Matlab of {variable: string, file: string}
			| Text of string
			| Dual of {text: string, matlab: {variable: string, file: string}}
		fun output genome =
			let
				val common = [
					("input", Misc.basenameWithoutExtension input)
					, ("order", Int.toString order)
					, ("genome", genome)
				]
				fun textName () =
					case Misc.substitute
						(("extension", "txt.gz") :: common) (!outFile)
					of
						NONE => Fail.fail "bad text output filename syntax"
						| SOME x => x
				fun matlabName () =
					case Misc.substitute
						(("extension", "mat") :: common) (!outFile)
					of
						NONE => Fail.fail "bad Matlab output filename syntax"
						| SOME x => x
			in
				case !outFormat of
					Otext => Text (textName ())
					| Omatlab => Matlab {
						variable = !variable
						, file = matlabName ()
					} | Odual => Dual {
						text = textName ()
						, matlab = {
							variable = !variable
							, file = matlabName ()
						}
					}
			end
	end
end
