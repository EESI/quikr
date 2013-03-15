signature PROGRAM = sig
	val version: string
	val name: string
end

structure Program :> PROGRAM = struct
	val version = "2.0"
	val name = "Naive Bayes Classifier - Score"
end
