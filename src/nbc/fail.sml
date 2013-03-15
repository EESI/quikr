signature FAIL = sig
	val fail: string -> 'a
end

structure Fail :> FAIL = struct
	fun fail why = (
		TextIO.output (TextIO.stdErr, why ^ "\n")
		; OS.Process.exit OS.Process.failure
	)
end
