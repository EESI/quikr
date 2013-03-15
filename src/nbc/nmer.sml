signature NMER_SIDES = sig
	type sidesBase
	type sidesNmer
	type sides
	val create: unit -> sides
	val clear: sides -> unit
	val put: sides * sidesBase -> unit
	val isFull: sides -> bool
	val forward: sides -> sidesNmer
end

signature NMER = sig
	eqtype base
	val a: base
	val c: base
	val g: base
	val t: base
	eqtype nmer
	val compare: nmer * nmer -> order
	val maximum: nmer
	val minimum: nmer
	val next: nmer -> nmer
	val hash: nmer -> Word.word
	val equal: nmer * nmer -> bool
	val toString: nmer -> string
	val fromString: string -> nmer option
	structure Single: NMER_SIDES
		where type sidesBase = base
		where type sidesNmer = nmer
	structure Double: sig
		include NMER_SIDES
			where type sidesBase = base
			where type sidesNmer = nmer
		val reverse: sides -> nmer
	end
end

signature NMER_ARGUMENTS = sig
	val order: int
	structure Word: sig
		eqtype word
		val fromInt: Int.int -> word
		val toInt: word -> Int.int
		val + : word * word -> word
		val << : word * Word.word -> word
		val ~>> : word * Word.word -> word
		val andb: word * word -> word
		val orb: word * word -> word
		val xorb: word * word -> word
		val compare: word * word -> order
		val toLarge: word -> LargeWord.word
	end
end

functor Nmer (Arguments: NMER_ARGUMENTS) = struct
	type base = Arguments.Word.word
	val a = Arguments.Word.fromInt 0
	val c = Arguments.Word.fromInt 1
	val g = Arguments.Word.fromInt 2
	val t = Arguments.Word.fromInt 3
	val maximumBase = t
	val baseBits = 0w2
	val nmerBits = Word.fromInt (Arguments.order * 2)
	fun opposite base =
		(*
			Conveniently enough, xor properly implements this:
				a -> t
				c -> g
				g -> c
				t -> a
		*)
		Arguments.Word.xorb (base, maximumBase)
	type nmer = Arguments.Word.word
	val compare = Arguments.Word.compare
	val minimum = Arguments.Word.fromInt 0
	local
		fun shiftInto (nmer, base) =
			Arguments.Word.+ (
				Arguments.Word.<< (nmer, baseBits)
				, base
			)
		fun maximumOfOrder order =
			if order = 0 then minimum
			else shiftInto (
				maximumOfOrder (order - 1)
				, maximumBase
			)
	in
		val maximum = maximumOfOrder Arguments.order
	end
	local
		val one = Arguments.Word.fromInt 1
	in
		fun next nmer = Arguments.Word.+ (nmer, one)
	end
	local
		fun charFromBase base = case Arguments.Word.toInt base of
			0 => #"A"
			| 1 => #"C"
			| 2 => #"G"
			| 3 => #"T"
			| _ => raise Fail "bug in nmer.sml"
		fun get (nmer, index) =
			let
				fun multiplyByTwo word = Word.<< (word, 0w1)
				val offset = multiplyByTwo (
					Word.fromInt (
						Arguments.order - 1 - index
					)
				)
			in
				Arguments.Word.~>> (
					Arguments.Word.andb (
						nmer
						, Arguments.Word.<< (
							maximumBase
							, offset
						)
					), offset
				)
			end
	in
		fun toString nmer = CharVector.tabulate (
			Arguments.order
			, fn index => charFromBase (
				get (nmer, index)
			)
		)
	end
	fun hash nmer = Word.fromLarge (Arguments.Word.toLarge nmer)
	fun equal (a, b) = Arguments.Word.compare (a, b) = EQUAL
	structure Undetermined = struct
		type sidesBase = base
		type sidesNmer = nmer
		type 'reverse undeterminedSides = {
			forward: nmer ref
			, reverse: 'reverse
			, count: int ref
		}
		fun clear {forward = _, reverse = _, count} = count := 0
		fun put ({forward, reverse, count}, base) = (
			forward := Arguments.Word.+ (
				Arguments.Word.andb (
					Arguments.Word.<< (
						!forward
						, baseBits
					), maximum
				), base
			);
				if !count = Arguments.order then ()
				else count := !count + 1
		)
		fun isFull {forward = _, reverse = _, count = ref count} =
			count = Arguments.order
		fun forward {forward = ref forward, reverse = _, count = _} =
			forward
	end
	structure Single = struct
		open Undetermined
		type sides = unit undeterminedSides
		fun create () = {
			forward = ref minimum
			, reverse = ()
			, count = ref 0
		}
	end
	structure Double = struct
		open Undetermined
		type sides = nmer ref undeterminedSides
		fun create () = {
			forward = ref minimum
			, reverse = ref maximum
			, count = ref 0
		}
		val put = fn (
			sides as {forward = _, reverse, count = _}
			, base
		) => (
			put (sides, base)
			; reverse := Arguments.Word.+ (
				Arguments.Word.~>> (
					!reverse
					, baseBits
				), Arguments.Word.<< (
					opposite base
					, nmerBits - baseBits
				)
			)
		)
		fun reverse {reverse = ref reverse, forward = _, count = _} =
			reverse
	end
	fun fromString string =
		let
			val side = Single.create ()
			val char = fn 
				#"A" => (Single.put (side, a); true)
				| #"a" => (Single.put (side, a); true)
				| #"C" => (Single.put (side, c); true)
				| #"c" => (Single.put (side, c); true)
				| #"G" => (Single.put (side, g); true)
				| #"g" => (Single.put (side, g); true)
				| #"T" => (Single.put (side, t); true)
				| #"t" => (Single.put (side, t); true)
				| _ => false
				
		in
			if CharVector.all char string then
				SOME (Single.forward side)
			else NONE
		end
end

functor Test () = struct
	structure Nmer1 = Nmer (
		val order = 1
		structure Word = Word32
	)
	structure Nmer2 = Nmer (
		val order = 2
		structure Word = Word32
	)
	structure Nmer3 = Nmer (
		val order = 3
		structure Word = Word32
	)
	val () = Test.list [
		{
			description = "opposite a = t"
			, function = fn () =>
				Bool.toString (
					Nmer1.opposite Nmer1.a = Nmer1.t
				)
			, expectedResult = "true"
		}, {
			description = "opposite c = g"
			, function = fn () =>
				Bool.toString (
					Nmer1.opposite Nmer1.c = Nmer1.g
				)
			, expectedResult = "true"
		}, {
			description = "A forward"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.a)
					; Nmer1.toString (
						Nmer1.Double.forward nmer
					)
				end
			, expectedResult = "A"
		}, {
			description = "A reverse"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.a)
					; Nmer1.toString (
						Nmer1.Double.reverse nmer
					)
				end
			, expectedResult = "T"
		}, {
			description = "C forward"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.c)
					; Nmer1.toString (
						Nmer1.Double.forward nmer
					)
				end
			, expectedResult = "C"
		}, {
			description = "C reverse"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.c)
					; Nmer1.toString (
						Nmer1.Double.reverse nmer
					)
				end
			, expectedResult = "G"
		}, {
			description = "G forward"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.g)
					; Nmer1.toString (
						Nmer1.Double.forward nmer
					)
				end
			, expectedResult = "G"
		}, {
			description = "G reverse"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.g)
					; Nmer1.toString (
						Nmer1.Double.reverse nmer
					)
				end
			, expectedResult = "C"
		}, {
			description = "T forward"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.t)
					; Nmer1.toString (
						Nmer1.Double.forward nmer
					)
				end
			, expectedResult = "T"
		}, {
			description = "T reverse"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.t)
					; Nmer1.toString (
						Nmer1.Double.reverse nmer
					)
				end
			, expectedResult = "A"
		}, {
			description = "AA forward"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.toString (
						Nmer2.Double.forward nmer
					)
				end
			, expectedResult = "AA"
		}, {
			description = "AA reverse"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.toString (
						Nmer2.Double.reverse nmer
					)
				end
			, expectedResult = "TT"
		}, {
			description = "AC forward"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.toString (
						Nmer2.Double.forward nmer
					)
				end
			, expectedResult = "AC"
		}, {
			description = "AC reverse"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.toString (
						Nmer2.Double.reverse nmer
					)
				end
			, expectedResult = "GT"
		}, {
			description = "GTA forward"
			, function = fn () =>
				let
					val nmer = Nmer3.Double.create ()
				in
					Nmer3.Double.put (nmer, Nmer3.g)
					; Nmer3.Double.put (nmer, Nmer3.t)
					; Nmer3.Double.put (nmer, Nmer3.a)
					; Nmer3.toString (
						Nmer3.Double.forward nmer
					)
				end
			, expectedResult = "GTA"
		}, {
			description = "GTA reverse"
			, function = fn () =>
				let
					val nmer = Nmer3.Double.create ()
				in
					Nmer3.Double.put (nmer, Nmer3.g)
					; Nmer3.Double.put (nmer, Nmer3.t)
					; Nmer3.Double.put (nmer, Nmer3.a)
					; Nmer3.toString (
						Nmer3.Double.reverse nmer
					)
				end
			, expectedResult = "TAC"
		}, {
			description = "( ) isFull"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Bool.toString (
						Nmer1.Double.isFull nmer
					)
				end
			, expectedResult = "false"
		}, {
			description = "(C) isFull"
			, function = fn () =>
				let
					val nmer = Nmer1.Double.create ()
				in
					Nmer1.Double.put (nmer, Nmer1.g)
					; Bool.toString (
						Nmer1.Double.isFull nmer
					)
				end
			, expectedResult = "true"
		}, {
			description = "(C ) isFull"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.c)
					; Bool.toString (
						Nmer2.Double.isFull nmer
					)
				end
			, expectedResult = "false"
		}, {
			description = "(CG) isFull"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.Double.put (nmer, Nmer2.g)
					; Bool.toString (
						Nmer2.Double.isFull nmer
					)
				end
			, expectedResult = "true"
		}, {
			description = "C(GA) isFull"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.Double.put (nmer, Nmer2.g)
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Bool.toString (
						Nmer2.Double.isFull nmer
					)
				end
			, expectedResult = "true"
		}, {
			description = "CGA(  ) isFull"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.Double.put (nmer, Nmer2.g)
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.clear nmer
					; Bool.toString (
						Nmer2.Double.isFull nmer
					)
				end
			, expectedResult = "false"
		}, {
			description = "CGA (AC) isFull"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.Double.put (nmer, Nmer2.g)
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.clear nmer
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.put (nmer, Nmer2.c)
					; Bool.toString (
						Nmer2.Double.isFull nmer
					)
				end
			, expectedResult = "true"
		}, {
			description = "CGA (AC) forward"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.Double.put (nmer, Nmer2.g)
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.clear nmer
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.toString (
						Nmer2.Double.forward nmer
					)
				end
			, expectedResult = "AC"
		}, {
			description = "CGA (AC) reverse"
			, function = fn () =>
				let
					val nmer = Nmer2.Double.create ()
				in
					Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.Double.put (nmer, Nmer2.g)
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.clear nmer
					; Nmer2.Double.put (nmer, Nmer2.a)
					; Nmer2.Double.put (nmer, Nmer2.c)
					; Nmer2.toString (
						Nmer2.Double.reverse nmer
					)
				end
			, expectedResult = "GT"
		}, {
			description = "TG fromString"
			, function = fn () =>
				case Nmer2.fromString "TG" of
					NONE => "invalid string"
					| SOME nmer => Nmer2.toString nmer
			, expectedResult = "TG"
		}
	]
end

functor Nmer (Arguments: NMER_ARGUMENTS) :> NMER = Nmer (Arguments)
