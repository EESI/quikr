signature BINARY = sig
	val fromInt32: int -> Word8Vector.vector
	val fromInt16: int -> Word8Vector.vector
	val fromReal: real -> Word8Vector.vector
end

structure Binary :> BINARY = struct
	val word8VectorFromArray = Word8ArraySlice.vector o Word8ArraySlice.full
	fun fromInt32 i =
		let
			val array = Word8Array.array (PackWord32Little.bytesPerElem, 0w0)
			val word = LargeWord.fromInt i
		in
			PackWord32Little.update (array, 0, word)
			; word8VectorFromArray array
		end
	fun fromInt16 i =
		let
			val array = Word8Array.array (PackWord16Little.bytesPerElem, 0w0)
			val word = LargeWord.fromInt i
		in
			PackWord16Little.update (array, 0, word)
			; word8VectorFromArray array
		end
	val fromReal = PackRealLittle.toBytes
end
