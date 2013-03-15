signature SUBSTITUTION = sig
	val substitute: (string -> string) -> string -> string option
end

structure Substitution :> SUBSTITUTION = struct
	structure LrVals = SubstitutionGrmLrValsFun (structure Token = LrParser.Token)
	structure Lex = SubstitutionLexFun (structure Tokens = LrVals.Tokens)
	structure Parser = JoinWithArg (
		structure ParserData = LrVals.ParserData
		structure Lex = Lex
		structure LrParser = LrParser
	)
	fun substitute lookup source =
		let
			val position = ref 0
			val instream = TextIO.openString source
			fun read n = TextIO.inputN (instream, n)
			val lexer = Parser.makeLexer read (ref 0)
			fun error (_, _, _) = ()
			val (result, _) = Parser.parse (0, lexer, error, lookup)
			val () = TextIO.closeIn instream
		in
			SOME result
		end
		handle Parser.ParseError => NONE
end
