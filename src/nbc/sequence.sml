signature SEQUENCE = sig
	type 'a t
	val fold: ('a * 'b -> 'b) -> 'b -> 'a t -> 'b
	val from: (unit -> 'a option) -> 'a t
	val map: ('a -> 'b) -> 'a t -> 'b t
	val app: ('a -> unit) -> 'a t -> unit
	val fromArray: 'a array -> 'a t
	val fromList: 'a list -> 'a t
	val toList: 'a t -> 'a list
	val toArray: 'a t -> 'a array
end

structure Sequence :> SEQUENCE = struct
	type 'a t = unit -> 'a option
	fun fold f seed t =
		let
			fun loop x =
				case t () of
					NONE => x
					| SOME y => loop (f (y, x))
		in
			loop seed
		end
	fun from t = t
	fun map f t =
		fn () => (
			case t () of
				NONE => NONE
				| SOME x => SOME (f x)
		)
	fun app f t =
		let
			fun loop () =
				case t () of
					NONE => ()
					| SOME x => (
						f x
						; loop ()
					)
		in
			loop ()
		end
	fun fromArray a =
		let
			val i = ref 0
			fun f () =
				if !i >= Array.length a then NONE
				else
					SOME (Array.sub (a, !i))
					before i := !i + 1
		in
			from f
		end
	fun fromList l =
		let
			val c = ref l
			fun f () =
				case !c of
					x :: y => (c := y; SOME x)
					| nil => NONE
		in
			from f
		end
	fun toList t =
		let
			fun loop l =
				case t () of
					NONE => rev l
					| SOME x => loop (x :: l)
		in
			loop nil
		end
	fun toArray t = Array.fromList (toList t)
end
