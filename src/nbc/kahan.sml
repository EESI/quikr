(* Kahan summation *)

signature KAHAN = sig
	type t
	val zero: t
	val add: t * real -> t
	val sum: t -> real
	val sequence: real Sequence.t -> real
	val list: real list -> real
	val array: real array -> real
end

structure Kahan :> KAHAN = struct
	type t = real * real
	val zero = (0.0, 0.0)
	fun add ((s, c), x) =
		let
			val y = x - c
			val t = s + y
		in
			(t, t - s - y)
		end
	fun sum (s, c) = s
	local
		fun swappedAdd (a, b) = add (b, a)
	in
		fun sequence e = sum (Sequence.fold swappedAdd zero e)
		fun list l = sum (foldl swappedAdd zero l)
		fun array a = sum (Array.foldl swappedAdd zero a)
	end
end
