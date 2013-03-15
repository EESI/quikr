signature JUDY = sig
	exception OutOfMemory
	type t
	val create: unit -> t
	val insert: t * string * int -> unit
	val get: t * string -> int option
	val bump: t * string -> unit
	val first: t -> (string * int) option
	val next: t * string -> (string * int) option
	val sequence: t -> (string * int) Sequence.t
	val app: (string * int -> unit) -> t -> unit
end
structure Judy :> JUDY = struct
	structure Primitive :> sig
		type judy
		type errorDetail
		type return
		val judyNull: judy
		val errorDetailNull: errorDetail
		val get: judy * string * errorDetail -> return
		val insert: judy ref * string * errorDetail -> return
		val delete: judy ref * string * errorDetail -> int
		val free: judy ref * errorDetail -> word
		val first: judy * CharArray.array * errorDetail -> return
		val next: judy * CharArray.array * errorDetail -> return
		val returnIsError: return -> bool
		val returnIsNull: return -> bool
		val returnGet: return -> int
		val returnSet: return * int -> unit
	end = struct
		type judy = MLton.Pointer.t
		type errorDetail = MLton.Pointer.t
		type return = MLton.Pointer.t
		val judyNull = MLton.Pointer.null
		val errorDetailNull = MLton.Pointer.null
		val get = _import "JudySLGet": judy * string * errorDetail -> return;
		val insert = _import "JudySLIns": judy ref * string * errorDetail -> return;
		val delete = _import "JudySLDel": judy ref * string * errorDetail -> int;
		val free = _import "JudySLFreeArray": judy ref * errorDetail -> word;
		val first = _import "JudySLFirst": judy * CharArray.array * errorDetail -> return;
		val next = _import "JudySLNext": judy * CharArray.array * errorDetail -> return;
		local
			val pjerr = MLton.Pointer.sub (MLton.Pointer.null, 0w1)
		in
			fun returnIsError return = return = pjerr
		end
		fun returnIsNull return = return = MLton.Pointer.null
		fun returnGet return = Int32.toInt (MLton.Pointer.getInt32 (return, 0))
		fun returnSet (return, i) = MLton.Pointer.setInt32 (return, 0, Int32.fromInt i)
	end
	exception OutOfMemory
	type t = {judy: Primitive.judy ref, max: int ref}
	fun create () = {judy = ref Primitive.judyNull, max = ref 0}
	fun insert ({judy, max}, key, value) =
		let
			val return = Primitive.insert (
				judy
				, key ^ "\000"
				, Primitive.errorDetailNull
			)
		in
			if Primitive.returnIsError return then raise OutOfMemory
			else let
				val n = size key
			in
				if !max < n then max := n else ()
				; Primitive.returnSet (return, value)
			end
		end
	fun get ({judy, max = _}, key) =
		let
			val return = Primitive.get (
				!judy
				, key ^ "\000"
				, Primitive.errorDetailNull
			)
		in
			if Primitive.returnIsNull return then NONE
			else SOME (Primitive.returnGet return)
		end
	fun bump ({judy, max}, key) =
		let
			val return = Primitive.insert (
				judy
				, key ^ "\000"
				, Primitive.errorDetailNull
			)
		in
			if Primitive.returnIsError return then raise OutOfMemory
			else let
				val n = size key
			in
				if !max < n then max := n else ()
				; Primitive.returnSet (
					return
					, Primitive.returnGet return + 1
				)
			end
		end
	fun strlen array =
		case CharArray.findi (fn (_, c) => c = #"\000") array of
			NONE => raise Option
			| SOME (i, _) => i
	fun stringFromNullTerminatedArray array =
		CharArraySlice.vector (
			CharArraySlice.slice (array, 0, SOME (strlen array))
		)
	fun first {judy, max} =
		let
			val array = CharArray.array (!max + 1, #"\000")
			val return = Primitive.first (
				!judy
				, array
				, Primitive.errorDetailNull
			)
		in
			if Primitive.returnIsNull return then NONE
			else SOME (
				stringFromNullTerminatedArray array
				, Primitive.returnGet return
			)
		end
	fun next ({judy, max}, key) =
		let
			val size = size key
			val array = CharArray.tabulate (
				!max + 1
				, fn i =>
					if i < size then String.sub (key, i)
					else #"\000"
			)
			val return = Primitive.next (
				!judy
				, array
				, Primitive.errorDetailNull
			)
		in
			if Primitive.returnIsNull return then NONE
			else SOME (
				stringFromNullTerminatedArray array
				, Primitive.returnGet return
			)
		end
	fun sequence t =
		let
			val last = ref NONE
			fun get () =
				case (
					case !last of
						NONE => first t
						| SOME key => next (t, key)
				) of
					NONE => NONE
					| SOME (return as (key, _)) => (
						last := SOME key
						; SOME return
					)
		in
			Sequence.from get
		end
	fun app f t =
		let
			fun apply (key, value) = (
				f (key, value)
				; fetch key
			) and fetch key =
				case next (t, key) of
					NONE => ()
					| SOME x => apply x
		in
			case first t of
				NONE => ()
				| SOME x => apply x
		end
end
