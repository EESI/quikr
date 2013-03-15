signature STREAM = sig
	type 'element stream
	val create:
		(unit -> ('element * 'element stream) option)
		-> 'element stream
	val empty: unit -> 'element stream
	val cons: 'element * 'element stream -> 'element stream
	val unfold:
		('seed -> ('fruit * 'seed) option)
		-> 'seed
		-> 'fruit stream
	val getItem: 'element stream -> ('element * 'element stream) option
	val isEmpty: 'element stream -> bool
	val fold:
		('element * 'accumulation -> 'accumulation)
		-> 'accumulation
		-> 'element stream
		-> 'accumulation
	val length: 'element stream -> int
	val rev: 'element stream -> 'element stream
	val map: ('input -> 'output) -> 'input stream -> 'output stream
	val mapPartial:
		('input -> 'output option)
		-> 'input stream
		-> 'output stream
	val app: ('element -> unit) -> 'element stream -> unit
	val toList: 'element stream -> 'element list
	val fromList: 'element list -> 'element stream
	val toVector: 'element stream -> 'element vector
	val fromVector: 'element vector -> 'element stream
	val fromVectorSlice: 'element VectorSlice.slice -> 'element stream
	val toArray: 'element stream -> 'element array
	val fromArray: 'element array -> 'element stream
	val fromArraySlice: 'element ArraySlice.slice -> 'element stream
	val fromString: string -> char stream
	val fromSubstring: Substring.substring -> char stream
	val toString: char stream -> string
	val fromTextInstream: TextIO.instream -> char stream
	val append: 'element stream * 'element stream -> 'element stream
	val concat: 'element stream stream -> 'element stream
	val hd: 'element stream -> 'element
	val tl: 'element stream -> 'element stream
	val find: ('element -> bool) -> 'element stream -> 'element option
	val filter: ('element -> bool) -> 'element stream -> 'element stream
	val exists: ('element -> bool) -> 'element stream -> bool
	val all: ('element -> bool) -> 'element stream -> bool
	val partition:
		('element -> bool)
		-> 'element stream
		-> 'element stream * 'element stream
	val take: ('element -> bool) -> 'element stream -> 'element stream
	val drop: ('element -> bool) -> 'element stream -> 'element stream
	val split:
		('element -> bool)
		-> 'element stream
		-> 'element stream * 'element stream
	val trim: 'element stream * int -> 'element stream
	val tokens:
		('element -> bool)
		-> 'element stream
		-> 'element stream stream
	val fields:
		('element -> bool)
		-> 'element stream
		-> 'element stream stream
end

structure Stream :> STREAM = struct
	datatype 'element stream =
		T of unit -> ('element * 'element stream) option
	fun create function = T function
	fun empty () = create (fn () => NONE)
	fun cons headAndTail = create (fn () => SOME headAndTail)
	fun unfold harvest seed = create (fn () =>
		case harvest seed of
			NONE => NONE
			| SOME (fruit, seed) => SOME (
				fruit
				, unfold harvest seed
			)
	)
	fun getItem (T function) = function ()
	fun fromList list = unfold List.getItem list
	fun toList stream = case getItem stream of
		NONE => nil
		| SOME (head, tail) => head :: toList tail
	fun fold accumulate accumulation stream = case getItem stream of
		NONE => accumulation
		| SOME (head, tail) =>
			fold accumulate (accumulate (head, accumulation)) tail
	fun length stream = fold (fn (_, count) => count + 1) 0 stream
	fun rev stream = fromList (fold op :: nil stream)
	fun map transform stream = unfold (fn stream => 
		case getItem stream of
			NONE => NONE
			| SOME (head, tail) => SOME (
				transform head
				, tail
			)
	) stream
	fun app execute stream =
		fold (fn (element, ()) => execute element) () stream
	fun fromVectorSlice slice = unfold VectorSlice.getItem slice
	fun fromVector vector = fromVectorSlice (VectorSlice.full vector)
	fun fromArraySlice slice = unfold ArraySlice.getItem slice
	fun fromArray array = fromArraySlice (ArraySlice.full array)
	fun fromSubstring substring = unfold Substring.getc substring
	fun fromString string = fromSubstring (Substring.full string)
	local
		fun withTabulate tabulate stream =
			let
				val position = ref stream
			in
				tabulate (
					length stream
					, fn _ => case getItem (!position) of
						NONE => raise Fail "Stream"
						| SOME (head, tail) => (
							position := tail
							; head
						)
				)
			end
	in
		fun toVector stream = withTabulate Vector.tabulate stream
		fun toArray stream = withTabulate Array.tabulate stream
		fun toString stream = withTabulate CharVector.tabulate stream
	end
	fun fromTextInstream instream =
		unfold TextIO.StreamIO.input1 (TextIO.getInstream instream)
	fun append (first, second) = create (fn () =>
		case getItem first of
			NONE => getItem second
			| SOME (head, tail) => SOME (
				head
				, append (tail, second)
			)
	)
	fun concat streams = create (fn () =>
		case getItem streams of
			NONE => NONE
			| SOME (head, tail) =>
				getItem (append (head, concat tail))
	)
	fun hd stream = case getItem stream of
		NONE => raise Empty
		| SOME (head, _) => head
	fun tl stream = case getItem stream of
		NONE => raise Empty
		| SOME (_, tail) => tail
	fun last stream = hd (rev stream)
	fun drop (stream, count) =
		if count < 0 then raise Subscript
		else if count = 0 then stream
		else case getItem stream of
			NONE => raise Subscript
			| SOME (_, tail) => drop (tail, count - 1)
	fun nth streamAndOffset = case getItem (drop streamAndOffset) of
		NONE => raise Subscript
		| SOME (head, _) => head
	fun mapPartial transform stream = create (fn () =>
		case getItem stream of
			NONE => NONE
			| SOME (head, tail) => case transform head of
				NONE => getItem (mapPartial transform tail)
				| SOME element => SOME (
					element
					, mapPartial transform tail
				)
	)
	fun find test stream = case getItem stream of
		NONE => NONE
		| SOME (head, tail) =>
			if test head then SOME head
			else find test tail
	fun filter test stream = unfold (fn stream =>
		case getItem stream of
			NONE => NONE
			| someHeadAndTail as (SOME (head, tail)) =>
				if test head then someHeadAndTail
				else getItem (filter test tail)
	) stream
	fun exists test stream = case find test stream of
		NONE => false
		| SOME _ => true
	fun all test stream = not (exists (not o test) stream)
	fun partition test stream =
		let
			val withResult = map (fn element =>
				(test element, element)
			) stream
		in (
			mapPartial (fn (result, element) =>
				if result then SOME element
				else NONE
			) withResult
			, mapPartial (fn (result, element) =>
				if result then NONE
				else SOME element
			) withResult
		) end
	fun take test stream = create (fn () =>
		case getItem stream of
			NONE => NONE
			| SOME (head, tail) =>
				if test head then SOME (head, take test tail)
				else NONE
	)
	fun drop test stream = create (fn () =>
		case getItem stream of
			NONE => NONE
			| someHeadAndTail as (SOME (head, tail)) =>
				if test head then getItem (drop test tail)
				else someHeadAndTail
	)
	fun split test stream = (take test stream, drop test stream)
	fun trim (stream, count) =
		if count <= 0 then stream
		else create (fn () =>
			case getItem stream of
				NONE => NONE
				| SOME (_, tail) =>
					getItem (trim (tail, count - 1))
		)
	fun isEmpty stream = case getItem stream of
		NONE => true
		| SOME _ => false
	fun tokens isSeparator stream = unfold (fn stream =>
		let
			val skipped = drop isSeparator stream
		in
			if isEmpty skipped then NONE
			else SOME (split (not o isSeparator) skipped)
		end
	) stream
	fun fields isSeparator stream = unfold (fn stream =>
		if isEmpty stream then NONE
		else SOME (
			take (not o isSeparator) stream
			, trim (drop (not o isSeparator) stream, 1)
		)
	) stream
end
