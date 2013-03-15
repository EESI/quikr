datatype sides = Single | Double
datatype labeled = Labeled | Unlabeled
local
	(*
	val perWord = ref NONE
	val total = ref NONE
	*)
	val order = ref (SOME 15)
	val sides = ref Double
	val labeled = ref Labeled
	val optionsWithoutHelp = [
		(* {
			short = "w", long = ["per-word"]
			, desc = GetOpt.ReqArg (
				fn file => perWord := SOME file
				, "file"
			), help = "file to store per-word counts in"
		}, {
			short = "t", long = ["total"]
			, desc = GetOpt.ReqArg (
				fn file => total := SOME file
				, "file"
			), help = "file to store total count in"
		}, *) {
			short = "r", long = ["order"]
			, desc = GetOpt.ReqArg (
				fn size => order := Int.fromString size
				, "size"
			), help = "word size"
		}, {
			short = "1", long = ["single"]
			, desc = GetOpt.NoArg (fn () => sides := Single)
			, help = "only count one side"
		}, {
			short = "u", long = ["unlabeled"]
			, desc = GetOpt.NoArg (fn () => labeled := Unlabeled)
			, help = "emit counts for every possible nmer, without labels"
		}
	]
	fun usageString () = GetOpt.usageInfo {
		header = CommandLine.name () ^ " <options> <input FASTA file> ..."
		, options = optionsWithoutHelp
	} ^ "\n"
	datatype status = Success | Failure
	fun displayHelpAndExit status = (
		TextIO.output (
			TextIO.stdErr
			, usageString ()
		); OS.Process.exit (case status of
			Success => OS.Process.success
			| Failure => OS.Process.failure
		)
	)
	val options = {
		short = "h", long = ["help"]
		, desc = GetOpt.NoArg (fn () => displayHelpAndExit Success)
		, help = "display help"
	} :: optionsWithoutHelp
in
	val (_, files) = GetOpt.getOpt {
		argOrder = GetOpt.Permute
		, options = options
		, errFn = fn errorMessage => (
			TextIO.output (TextIO.stdErr, errorMessage ^ "\n")
			; displayHelpAndExit Failure
		)
	} (CommandLine.arguments ())
	(*
	val perWordFileName = case !perWord of
		NONE => (
			TextIO.output (
				stdErr
				, "per-word file name required but not provided\n"
			); displayHelpAndExit Failure
		) | SOME fileName => fileName
	val totalFileName = case !total of
		NONE => (
			TextIO.output (
				stdErr
				, "total file name required but not provided\n"
			); displayHelpAndExit Failure
		) | SOME fileName => fileName
	*)
	val order = case !order of
		NONE => (
			TextIO.output (
				TextIO.stdErr
				, "invalid order\n"
			); displayHelpAndExit Failure
		) | SOME integer => integer
	val sides = !sides
	val labeled = !labeled
end

signature COLLECTION = sig
	type collection
	type nmer
	val empty: unit -> collection
	val add: collection * nmer -> unit
	val get: collection * nmer -> int
	val app: (nmer * int -> unit) -> collection -> unit
end

functor Collection (Nmer: NMER)
:> COLLECTION where type nmer = Nmer.nmer = struct
	type nmer = Nmer.nmer
	structure Table = HashTableFn (
		type hash_key = nmer
		val hashVal = Nmer.hash
		val sameKey = Nmer.equal
	)
	type collection = int ref Table.hash_table
	exception NotFound
	fun empty () = Table.mkTable (256 * 1024, NotFound)
	fun add (table, nmer) = case Table.find table nmer of
		NONE => Table.insert table (nmer, ref 1)
		| SOME count => count := !count + 1
	fun get (table, nmer) = case Table.find table nmer of
		NONE => 0
		| SOME (ref count) => count
	fun app execute table = Table.appi (fn (nmer, ref count) =>
		execute (nmer, count)
	) table
end

datatype result = Success | Failure

signature OUTPUT = sig
	type collection
	val output: collection -> unit
end

functor Unlabeled (
	structure Nmer: NMER
	structure Collection: COLLECTION
		sharing type Collection.nmer = Nmer.nmer
) :> OUTPUT
	where type collection = Collection.collection
= struct
	type collection = Collection.collection
	fun put string = TextIO.output (TextIO.stdOut, string)
	fun single count = (
		put (Int.toString count)
		; put "\n"
	)
	fun output collection =
		let
			fun continue nmer = (
				single (Collection.get (collection, nmer))
				;
					if nmer = Nmer.maximum then ()
					else continue (Nmer.next nmer)
			)
		in
			continue (Nmer.minimum)
		end
end

functor Labeled (
	structure Nmer: NMER
	structure Collection: COLLECTION
		sharing type Collection.nmer = Nmer.nmer
) :> OUTPUT
	where type collection = Collection.collection
= struct
	type collection = Collection.collection
	fun put string = TextIO.output (TextIO.stdOut, string)
	fun single (nmer, count) = (
		put (Nmer.toString nmer)
		; put " "
		; put (Int.toString count)
		; put "\n"
	)
	fun output collection = Collection.app single collection
end

functor File (
	structure Collection: COLLECTION
	structure Output: OUTPUT
		sharing type Collection.collection = Output.collection
) :> FILE
	where type nmer = Collection.nmer
	where type result = result
	where type argument = unit
= struct
	type argument = unit
	type file = Collection.collection
	type read = unit
	type nmer = Collection.nmer
	type result = result
	fun startFile _ = Collection.empty ()
	fun startRead _ = ()
	fun nmer (counts, (), nmer) = Collection.add (counts, nmer)
	fun stopRead (_, ()) = ()
	fun stopFile counts = (
		Output.output counts
		; Success
	)
	fun invalidFormat file = Failure
end

functor Everything (Nmer: NMER) = struct
	structure Collection = Collection (Nmer)
	structure Unlabeled = File (
		structure Collection = Collection
		structure Output = Unlabeled (
			structure Nmer = Nmer
			structure Collection = Collection
		)
	)
	structure Labeled = File (
		structure Collection = Collection
		structure Output = Labeled (
			structure Nmer = Nmer
			structure Collection = Collection
		)
	)
	structure SingleSidedUnlabeled = SingleSidedFasta (
		structure Nmer = Nmer
		structure File = Unlabeled
	)
	structure DoubleSidedUnlabeled = DoubleSidedFasta (
		structure Nmer = Nmer
		structure File = Unlabeled
	)
	structure SingleSidedLabeled = SingleSidedFasta (
		structure Nmer = Nmer
		structure File = Labeled
	)
	structure DoubleSidedLabeled = DoubleSidedFasta (
		structure Nmer = Nmer
		structure File = Labeled
	)
end

structure Everything32 = Everything (
	Nmer (
		val order = order
		structure Word = Word32
	)
)
structure Everything64 = Everything (
	Nmer (
		val order = order
		structure Word = Word64
	)
)

val process =
	if order <= 32 then (case sides of
		Single => (case labeled of
			Unlabeled => Everything32.SingleSidedUnlabeled.process
			| Labeled => Everything32.SingleSidedLabeled.process
		) | Double => (case labeled of
			Unlabeled => Everything32.DoubleSidedUnlabeled.process
			| Labeled => Everything32.DoubleSidedLabeled.process
		)
	) else (case sides of
		Single => (case labeled of
			Unlabeled => Everything64.SingleSidedUnlabeled.process
			| Labeled => Everything64.SingleSidedLabeled.process
		) | Double => (case labeled of
			Unlabeled => Everything64.DoubleSidedUnlabeled.process
			| Labeled => Everything64.DoubleSidedLabeled.process
		)
	)

val () =
	let
		fun one name =
			let
				val instream = TextIO.openIn name
				val result = process ((), instream)
			in
				TextIO.closeIn instream
				; case result of
					Success => true
					| Failure => (
						TextIO.output (
							TextIO.stdErr
							, name
							^ ": invalid format\n"
						); false
					)
			end
		fun all names = List.all one names
	in
		if all files then ()
		else OS.Process.exit OS.Process.failure
	end
