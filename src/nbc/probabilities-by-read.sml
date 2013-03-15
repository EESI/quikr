local
	fun usage result = (
		TextIO.output (
			TextIO.stdErr
			, CommandLine.name ()
		^ " <order> <input FASTA> <file of nmers to count>\n"
		); OS.Process.exit result
	)
in
	val (order, input, toCount) = case CommandLine.arguments () of
		[order, input, toCount] => (case Int.fromString order of
			NONE => usage OS.Process.failure
			| SOME order => (order, input, toCount)
		) | ["--help"] => usage OS.Process.success
		| _ => usage OS.Process.failure
end

datatype result = Success | Failure
fun warn message = TextIO.output (
	TextIO.stdErr
	, (
		"warning: "
		^ message
		^ "\n"
	)
)

signature NMER_TABLE = sig
	type nmer
	type table

	(*
		Create the table. Only the provided nmers will be counted.
	*)
	val create: nmer list -> table

	(*
		Increment the count for a given nmer. If the nmer was
		not provided at table creation, do nothing.
	*)
	val bump: table * nmer -> unit

	(*
		Reset all counts to zero. Do not change the list of
		nmers to count.
	*)
	val clear: table -> unit

	(*
		Apply a function to all nmers and their counts, in
		lexicographic order.
	*)
	val app: (nmer * int -> unit) -> table -> unit
end

functor NmerTable (Nmer: NMER)
:> NMER_TABLE where type nmer = Nmer.nmer
= struct
	structure HashTable = HashTableFn (
		type hash_key = Nmer.nmer
		val hashVal = Nmer.hash
		val sameKey = Nmer.equal
	)
	exception NotFound
	type nmer = Nmer.nmer
	type table = {
		indexes: int HashTable.hash_table
		, counts: int array
		, nmers: nmer vector
	}
	fun create list =
		let
			val indexes = HashTable.mkTable (1024, NotFound)
			val nmers =
				let
					val array = Array.fromList list
				in
					ArrayQSort.sort Nmer.compare array
					; Array.vector array
				end
		in
			Vector.appi (fn (index, nmer) =>
				HashTable.insert indexes (nmer, index)
			) nmers
			; {
				indexes = indexes
				, nmers = nmers
				, counts = Array.array (
					Vector.length nmers
					, 0
				)
			}
		end
	fun bump ({indexes, nmers = _, counts}, nmer) =
		case HashTable.find indexes nmer of
			NONE => ()
			| SOME index => Array.update (
				counts
				, index
				, Array.sub (counts, index) + 1
			)
	fun clear {indexes = _, nmers = _, counts} =
		Array.modify (fn _ => 0) counts
	fun app execute {indexes = _, nmers, counts} =
		Vector.appi (fn (index, nmer) =>
			execute (nmer, Array.sub (counts, index))
		) nmers
end

functor Input (
	structure Nmer: NMER
	structure NmerTable: NMER_TABLE
	sharing type Nmer.nmer = NmerTable.nmer
) = SingleSidedFasta (
	structure Nmer = Nmer
	structure File = struct
		type argument = NmerTable.table
		type file = argument
		type read = Int64.int ref
		type nmer = Nmer.nmer
		datatype result = datatype result
		exception NotFound
		fun startFile table = table
		fun startRead (table, _) = ref (0: Int64.int)
		fun nmer (table, (total: Int64.int ref), nmer) = (
			total := !total + 1
			; NmerTable.bump (table, nmer)
		)
		fun put string = TextIO.output (TextIO.stdOut, string)
		fun stopRead (table, total) =
			let
				val realFromInt64 =
					Real.fromLargeInt o Int64.toLarge
				val realTotal = realFromInt64 (!total)
				val toString = Real.fmt (
					StringCvt.FIX (SOME 17)
				)
				fun probability count = real count / realTotal
				infix |>
				fun argument |> function = function argument
				val first = ref true
			in
				NmerTable.app (fn (nmer, count) => (
					if !first then first := false
					else put "\t"
					; count
						|> probability
						|> toString
						|> put
				)) table
				; put "\n"
				; NmerTable.clear table
				; total := 0
			end
		fun stopFile _ = Success
		fun invalidFormat _ = Failure
	end
)

structure Nmer = Nmer (
	val order = order
	structure Word = Word64
)
structure NmerTable = NmerTable (Nmer)
structure Input = Input (
	structure Nmer = Nmer
	structure NmerTable = NmerTable
)

val table =
	let
		fun build collect goose = case collect goose of
			NONE => nil
			| SOME egg =>
				egg :: build collect goose
		fun chopEnd string = String.extract (
			string
			, 0
			, SOME (size string - (
				if String.isSuffix "\r\n" string
				then 2
				else 1
			))
		)
		val line = Option.map chopEnd o TextIO.inputLine
		val lines = build line
		val instream = TextIO.openIn toCount
		fun fromStringWithWarning string =
			case Nmer.fromString string of
				NONE => (
					warn (
						string
						^ " is not a valid nmer"
					); NONE
				) | someNmer => someNmer
		val nmers = List.mapPartial
			fromStringWithWarning
			(lines instream)
	in
		TextIO.closeIn instream
		; NmerTable.create nmers
	end

val () = case Input.process (table, TextIO.openIn input) of
	Failure => (
		TextIO.output (
			TextIO.stdErr
			, "input is not valid FASTA\n"
		); OS.Process.exit OS.Process.failure
	) | Success => ()
