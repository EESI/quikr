functor Tree (Key: sig
	type key
	val compare: key * key -> order
end) :> sig
	type key = Key.key
	type 'datum tree
	val size: 'datum tree -> int
	val empty: 'datum tree
	val single: key * 'datum -> 'datum tree
	val add: 'datum tree * key * 'datum -> 'datum tree
	val find: 'datum tree * key -> 'datum option
	val exists: 'datum tree * key -> bool
	val up: 'datum tree -> (key * 'datum) Stream.stream
	val down: 'datum tree -> (key * 'datum) Stream.stream
	val upFrom: 'datum tree * key -> (key * 'datum) Stream.stream
	val downFrom: 'datum tree * key -> (key * 'datum) Stream.stream
	val fromList: (key * 'datum) list -> 'datum tree
	val fromStream: (key * 'datum) Stream.stream -> 'datum tree
	(*
	val remove: 'datum tree * key -> 'datum tree
	*)
end = struct
	type key = Key.key
	structure Height = Int8
	datatype 'datum tree =
		Leaf
		| Branch of {
			height: Height.int
			, less: 'datum tree
			, key: key
			, datum: 'datum
			, greater: 'datum tree
		}
	fun size tree = case tree of
		Leaf => 0
		| Branch {less, greater, ...} =>
			1 + size less + size greater
	val empty = Leaf
	fun single (key, datum) = Branch {
		height = 1
		, less = Leaf
		, key = key
		, datum = datum
		, greater = Leaf
	}
	fun height tree = case tree of
		Leaf => 0
		| Branch {height, ...} => height
	fun calculateHeight {key, datum, less, greater} = Branch {
		key = key
		, datum = datum
		, less = less
		, greater = greater
		, height = 1 + Height.max (height less, height greater)
	}
	fun rotateLess branch = case branch of
	(*
		   b          d
		 a   d  =>  b   e
		    c e    a c
	*)
		{
			less = a
			, key = bKey
			, datum = bDatum
			, greater = Branch {
				less = c
				, key = dKey
				, datum = dDatum
				, greater = e
				, height = _
			}
		} => calculateHeight {
			less = calculateHeight {
				less = a
				, key = bKey
				, datum = bDatum
				, greater = c
			}, key = dKey
			, datum = dDatum
			, greater = e
		} | _ => raise Fail "rotateLess"
	fun rotateGreater branch = case branch of
	(*
		   d        b
		 b   e => a   d
		a c          c e
	*)
		{
			less = Branch {
				less = a
				, key = bKey
				, datum = bDatum
				, greater = c
				, height = _
			}, key = dKey
			, datum = dDatum
			, greater = e
		} => calculateHeight {
			less = a
			, key = bKey
			, datum = bDatum
			, greater = calculateHeight {
				less = c
				, key = dKey
				, datum = dDatum
				, greater = e
			}
		} | _ => raise Fail "rotateGreater"
	fun balance (branch as {key, datum, less, greater}) =
		let
			val heightLess = height less
			val heightGreater = height greater
		in
			if heightLess < heightGreater - 2 then
				rotateLess branch
			else if heightLess > heightGreater + 2 then
				rotateGreater branch
			else calculateHeight branch
		end
	fun add (tree, newKey, newDatum) = case tree of
		Leaf => single (newKey, newDatum)
		| Branch {height, less, key, datum, greater} =>
			case Key.compare (newKey, key) of
				EQUAL => Branch {
					height = height
					, less = less
					, key = newKey
					, datum = newDatum
					, greater = greater
				} | LESS => balance {
					less = add (less, newKey, newDatum)
					, key = key
					, datum = datum
					, greater = greater
				} | GREATER => balance {
					less = less
					, key = key
					, datum = datum
					, greater = add (greater, newKey, newDatum)
				}
	fun find (tree, desiredKey) = case tree of
		Leaf => NONE
		| Branch {less, key, datum, greater, ...} =>
			case Key.compare (desiredKey, key) of
				EQUAL => SOME datum
				| LESS => find (less, desiredKey)
				| GREATER => find (greater, desiredKey)
	fun exists (tree, desiredKey) = case find (tree, desiredKey) of
		NONE => false
		| SOME _ => true
	fun up tree = Stream.create (fn () =>
		case tree of
			Leaf => NONE
			| Branch {less, key, datum, greater, ...} =>
				Stream.getItem (
					Stream.append (
						up less
						, Stream.cons (
							(key, datum)
							, up greater
						)
					)
				)
	)
	fun down tree = Stream.create (fn () =>
		case tree of
			Leaf => NONE
			| Branch {greater, key, datum, less, ...} =>
				Stream.getItem (
					Stream.append (
						down greater
						, Stream.cons (
							(key, datum)
							, down less
						)
					)
				)
	)
	fun upFrom (tree, firstKey) = Stream.create (fn () =>
		case tree of
			Leaf => NONE
			| Branch {less, key, datum, greater, ...} =>
				case Key.compare (firstKey, key) of
					LESS => Stream.getItem (
						Stream.append (
							upFrom (less, firstKey)
							, Stream.cons (
								(key, datum)
								, up greater
							)
						)
					) | EQUAL => SOME (
						(key, datum)
						, up greater
					) | GREATER => Stream.getItem (
						upFrom (greater, firstKey)
					)
	)
	fun downFrom (tree, firstKey) = Stream.create (fn () =>
		case tree of
			Leaf => NONE
			| Branch {greater, key, datum, less, ...} =>
				case Key.compare (firstKey, key) of
					LESS => Stream.getItem (
						downFrom (less, firstKey)
					) | EQUAL => SOME (
						(key, datum)
						, down less
					) | GREATER => Stream.getItem (
						Stream.append (
							downFrom (greater, firstKey)
							, Stream.cons (
								(key, datum)
								, down less
							)
						)
					)
	)
	fun fromStream stream =
		Stream.fold (fn ((key, datum), tree) =>
			add (tree, key, datum)
		) empty stream
	fun fromList list = fromStream (Stream.fromList list)
end
