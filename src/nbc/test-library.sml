signature TEST = sig
	type test =
		{
			description: string
			, expectedResult: string
			, function: unit -> string
		}
	val single: test -> unit
	val list: test list -> unit
end

structure Test = struct
	fun single {description, expectedResult, function} =
		let
			val actualResult = function ()
		in
			if expectedResult = actualResult then
				TextIO.output (
					TextIO.stdErr
					, (
						description
						^ " succeeded.\n"
					)
				)
			else (
				TextIO.output (
					TextIO.stdErr
					, (
						description
						^ " was supposed to be "
						^ expectedResult
						^ ", but was actually "
						^ actualResult
						^ ".\n"
					)
				); OS.Process.exit OS.Process.failure
			)
		end handle exception' => (
			TextIO.output (
				TextIO.stdErr
				, (
					description
					^ " failed with exception "
					^ exnMessage exception'
					^ ".\n"
				)
			); OS.Process.exit OS.Process.failure
		)
	fun list tests = app single tests
end
