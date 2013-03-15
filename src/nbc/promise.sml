structure Promise
:> sig
	type 'fulfillment promise
	val delay: (unit -> 'fulfillment) -> 'fulfillment promise
	val force: 'fulfillment promise -> 'fulfillment
end = struct
	local
                datatype 'expectation lazy =
                        Delayed of unit -> 'expectation
                        | Forced of 'expectation
        in
                type 'expectation promise = 'expectation lazy ref
                fun delay fulfill = ref (Delayed fulfill)
                fun force promise = case !promise of
                        Delayed fulfill =>
                                let
                                        val expectation = fulfill ()
                                in
                                        promise := Forced expectation
                                        ; expectation
                                end
                        | Forced expectation => expectation
        end
end
