signature STORE_JUDY = sig
	type t
	val load: (int * string) Sequence.t -> t
	val get: t * string -> int option
end

structure StoreJudy :> STORE_JUDY = struct
	type t = Judy.t
	fun load e =
		let
			val j = Judy.create ()
		in
			Sequence.app (fn (count, nmer) => Judy.insert (j, nmer, count)) e
			; j
		end
	val get = Judy.get
end
