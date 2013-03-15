signature FILE = sig
	type argument
	type file
	type read
	type nmer
	type result
	val startFile: argument -> file
	val startRead: file * string -> read
	val nmer: file * read * nmer -> unit
	val stopRead: file * read -> unit
	val stopFile: file -> result
	val invalidFormat: file -> result
end

signature FASTA = sig
	type argument
	type result
	val process: argument * TextIO.instream -> result
end

functor AgnosticFasta (
	structure Nmer: NMER
	structure File: FILE
		sharing type Nmer.nmer = File.nmer
	structure Sides: sig
		include NMER_SIDES
		type file
		type read
		val process: file * read * sides -> unit
	end
		sharing type Nmer.base = Sides.sidesBase
		sharing type Nmer.nmer = Sides.sidesNmer
		sharing type File.read = Sides.read
		sharing type File.file = Sides.file
) :> FASTA
	where type argument = File.argument
	where type result = File.result
= struct
	type argument = File.argument
	type result = File.result

	val beforeHeaderBeginningOfLine = ParseState.create ()
	val beforeHeaderMiddleOfLine = ParseState.create ()
	val afterHeaderBeginningOfLine = ParseState.create ()
	val afterHeaderMiddleOfLine = ParseState.create ()

	fun inputLineButDiscardNewline instream =
		Option.map (fn line =>
			String.extract (line, 0, SOME (size line - 1))
		) (TextIO.inputLine instream)
	datatype z = datatype ParseState.whichCharacters (* This | Any *)

	local
		fun header (instream, (file, sides)) =
			case inputLineButDiscardNewline instream of
				NONE => File.invalidFormat file
				| SOME header => ParseState.enter (
					afterHeaderBeginningOfLine
					, instream
					, (
						file
						, File.startRead (
							file
							, header
						), sides
					)
				)
		fun space (instream, (file, sides)) = ParseState.enter (
			beforeHeaderMiddleOfLine
			, instream
			, (file, sides)
		)
		fun newline (instream, (file, sides)) = ParseState.enter (
			beforeHeaderBeginningOfLine
			, instream
			, (file, sides)
		)
		fun invalidFormat (_, (file, _)) = File.invalidFormat file
	in
		val () = ParseState.build {
			state = beforeHeaderBeginningOfLine
			, characters = [
				(These [#">"], header)
				, (These [#"\n"], newline)
				, (These [#" ", #"\t", #"\r"], space)
				, (Any, invalidFormat)
			], endOfFile = invalidFormat
		}
		val () = ParseState.build {
			state = beforeHeaderMiddleOfLine
			, characters = [
				(These [#"\n"], newline)
				, (These [#" ", #"\t", #"\r"], space)
				, (Any, invalidFormat)
			], endOfFile = invalidFormat
		}
	end
	local
		fun base base (instream, (file, read, sides)) = (
			Sides.put (sides, base)
			;
				if Sides.isFull sides then
					Sides.process (file, read, sides)
				else ()
			; ParseState.enter (
				afterHeaderMiddleOfLine
				, instream
				, (file, read, sides)
			)
		)
		fun space (instream, (file, read, sides)) = (
			ParseState.enter (
				afterHeaderMiddleOfLine
				, instream
				, (file, read, sides)
			)
		)
		fun other (instream, (file, read, sides)) = (
			Sides.clear sides
			; ParseState.enter (
				afterHeaderMiddleOfLine
				, instream
				, (file, read, sides)
			)
		)
		fun newline (instream, (file, read, sides)) =
			ParseState.enter (
				afterHeaderBeginningOfLine
				, instream
				, (file, read, sides)
			)
		fun header (instream, (file, read, sides)) = (
			File.stopRead (file, read)
			; Sides.clear sides
			; case inputLineButDiscardNewline instream of
				NONE => File.invalidFormat file
				| SOME header => ParseState.enter (
					afterHeaderBeginningOfLine
					, instream
					, (
						file
						, File.startRead (
							file
							, header
						), sides
					)
				)
		)
		fun success (_, (file, read, _)) = (
			File.stopRead (file, read)
			; File.stopFile file
		)
	in
		val () = ParseState.build {
			state = afterHeaderBeginningOfLine
			, characters = [
				(These [#"A", #"a"], base Nmer.a)
				, (These [#"C", #"c"], base Nmer.c)
				, (These [#"G", #"g"], base Nmer.g)
				, (These [#"T", #"t"], base Nmer.t)
				, (These [#">"], header)
				, (These [#"\n"], newline)
				, (These [#" ", #"\t", #"\r"], space)
				, (Any, other)
			], endOfFile = success
		}
		val () = ParseState.build {
			state = afterHeaderMiddleOfLine
			, characters = [
				(These [#"A", #"a"], base Nmer.a)
				, (These [#"C", #"c"], base Nmer.c)
				, (These [#"G", #"g"], base Nmer.g)
				, (These [#"T", #"t"], base Nmer.t)
				, (These [#" ", #"\t", #"\r"], space)
				, (These [#"\n"], newline)
				, (Any, other)
			], endOfFile = success
		}
	end
	fun process (argument, instream) = ParseState.enter (
		beforeHeaderBeginningOfLine
		, instream
		, (File.startFile argument, Sides.create ())
	)
end

functor SingleSidedFasta (
	structure Nmer: NMER
	structure File: FILE
		sharing type Nmer.nmer = File.nmer
) = AgnosticFasta (
	structure Nmer = Nmer
	structure File = File
	structure Sides = struct
		type read = File.read
		type file = File.file
		open Nmer.Single
		fun process (file, read, sides) =
			File.nmer (file, read, forward sides)
	end
)

functor DoubleSidedFasta (
	structure Nmer: NMER
	structure File: FILE
		sharing type Nmer.nmer = File.nmer
) = AgnosticFasta (
	structure Nmer = Nmer
	structure File = File
	structure Sides = struct
		type read = File.read
		type file = File.file
		open Nmer.Double
		fun process (file, read, sides) = (
			File.nmer (file, read, forward sides)
			; File.nmer (file, read, reverse sides)
		)
	end
)

functor TestFile (Nmer: NMER) = struct
	type argument = unit
	type nmer = Nmer.nmer
	type read = {header: string, nmers: nmer list ref}
	type file = {header: string, nmers: string list} list ref
	type result = string
	fun startFile () = ref nil
	fun stopFile file = String.concatWith ";" (
		map (fn {header, nmers} =>
			header
			^ ":"
			^ String.concatWith "," (rev nmers)
		) (rev (!file))
	)
	fun startRead (_, header) =
		{header = header, nmers = ref nil}
	fun nmer (_, {header = _, nmers}, nmer) =
		nmers := nmer :: !nmers
	fun stopRead (file, {header, nmers = ref nmers}) =
		file := {
			header = header
			, nmers = map Nmer.toString nmers
		} :: !file
	fun invalidFormat _ = "invalid format"
end

functor Test () = struct
	structure Nmer1 = Nmer (
		val order = 1
		structure Word = Word32
	)
	structure File1 = TestFile (Nmer1)
	structure SingleFasta1 = SingleSidedFasta (
		structure Nmer = Nmer1
		structure File = File1
	)
	fun test process input () = process ((), TextIO.openString input)
	val single1 = test SingleFasta1.process
	structure Nmer2 = Nmer (
		val order = 2
		structure Word = Word32
	)
	structure File2 = TestFile (Nmer2)
	structure SingleFasta2 = SingleSidedFasta (
		structure Nmer = Nmer2
		structure File = File2
	)
	val single2 = test SingleFasta2.process
	structure DoubleFasta1 = DoubleSidedFasta (
		structure Nmer = Nmer1
		structure File = File1
	)
	val double1 = test DoubleFasta1.process
	structure DoubleFasta2 = DoubleSidedFasta (
		structure Nmer = Nmer2
		structure File = File2
	)
	val double2 = test DoubleFasta2.process
	val () = Test.list [
		{
			description = "single 1: A"
			, function = single1 ">foo\nA\n"
			, expectedResult = "foo:A"
		}, {
			description = "single 1: AG"
			, function = single1 ">foo\nAG\n"
			, expectedResult = "foo:A,G"
		}, {
			description = "single 2: A"
			, function = single2 ">foo\nA\n"
			, expectedResult = "foo:"
		}, {
			description = "single 2: CTGAG"
			, function = single2 ">foo\nCTGAG\n"
			, expectedResult = "foo:CT,TG,GA,AG"
		}, {
			description = "double 1: C"
			, function = double1 ">bar\nC\n"
			, expectedResult = "bar:C,G"
		}, {
			description = "double 2: T"
			, function = double2 ">baz\nT\n"
			, expectedResult = "baz:"
		}, {
			description = "double 2: GC"
			, function = double2 ">quux\nGC\n"
			, expectedResult = "quux:GC,GC"
		}, {
			description = "double 2: CCC\\nC\\nCT"
			, function = double2 ">goo\nCCC\nC\nCT\n"
			, expectedResult = "goo:CC,GG,CC,GG,CC,GG,CC,GG,CT,AG"
		}, {
			description = "double 2: CC\\nC*\\nT"
			, function = double2 ">goo\nCC\nC*\nT\n"
			, expectedResult = "goo:CC,GG,CC,GG"
		}, {
			description = "double 2: foo CATGAC goo TACCAG"
			, function = double2
				">foo\nCATGAC\n>goo\nTACCAG\n"
			, expectedResult = (
				"foo:CA,TG,AT,AT,TG,CA,GA,TC,AC,GT"
				^ ";goo:TA,TA,AC,GT,CC,GG,CA,TG,AG,CT"
			)
		}
	]
end
