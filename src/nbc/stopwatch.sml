signature STOPWATCH = sig
	exception FinishWithoutStart
	val start: string -> unit
	val finish: unit -> unit
end
structure Stopwatch :> STOPWATCH = struct
	exception FinishWithoutStart
	local
		val time = ref NONE
	in
		fun start doing = (
			print (concat [doing, "..."])
			; time := SOME (Time.now ())
		)
		fun finish () = case !time of
			SOME t => (
				print (concat [
					" done in "
					, Time.toString (Time.- (Time.now (), t))
					, " seconds.\n"
				]); time := NONE
			) | NONE => raise FinishWithoutStart
	end
end
