signature MATLAB = sig
	type t
	val openOut: string -> t
	type doubleArray
	val beginDoubleArray: t * string -> doubleArray
	val writeDouble: doubleArray * real -> unit
	val concludeDoubleArray: doubleArray -> unit
	val writeDoubleArray: t * string * real vector -> unit
	val closeOut: t -> unit
end

structure Matlab :> MATLAB = struct
	fun outputReal (c, f) = BinIO.output (c, Binary.fromReal f)
	fun outputInt32 (c, i) = BinIO.output (c, Binary.fromInt32 i)
	fun outputInt16 (c, i) = BinIO.output (c, Binary.fromInt16 i)
	type t = BinIO.outstream
	type doubleArray = {
		outstream: BinIO.outstream
		, wholeSizePos: BinIO.StreamIO.out_pos
		, paddedSize: int
		, arraySizePos: BinIO.StreamIO.out_pos
		, dataSizePos: BinIO.StreamIO.out_pos
		, elements: int ref
	}
	fun pad n =
		if n > 0 andalso n <= 4 then 4
		else (n + 7) div 8 * 8
	fun doubleArraySize (nameLength, arrayLength) =
		32 + nameLength + 4 + 8 + arrayLength * 8
	val s8Tag = 1
	val s32Tag = 5
	val u32Tag = 6
	val doubleTag = 9
	val doubleFlag = 6
	val matrixTag = 14
	fun lsl (x, y) = Word.toIntX (Word.<< (Word.fromInt x, Word.fromInt y))
	fun beginDoubleArray (c, name) =
		let
			val nameSize = Int.min (size name, 63)
			val name = Word8VectorSlice.vector (
				Word8VectorSlice.slice (Byte.stringToBytes name, 0, SOME nameSize)
			)
			val paddedSize = pad nameSize
			val padding = Word8Vector.tabulate (paddedSize - nameSize, fn _ => 0w0)
			val () = outputInt32 (c, matrixTag)
			val wholeSizePos = BinIO.getPosOut c
			val () = (
				outputInt32 (c, 0)
				; outputInt32 (c, u32Tag)
				; outputInt32 (c, 8)
				; outputInt32 (c, doubleFlag)
				; outputInt32 (c, 0)
				; outputInt32 (c, s32Tag)
				; outputInt32 (c, 8)
				; outputInt32 (c, 1)
			)
			val arraySizePos = BinIO.getPosOut c
			val () = (
				outputInt32 (c, 0)
				; if nameSize > 0 andalso nameSize <= 4 then
					outputInt32 (c, lsl (nameSize, 16) + s8Tag)
				else (
					outputInt32 (c, s8Tag)
					; outputInt32 (c, nameSize)
				); BinIO.output (c, name)
				; BinIO.output (c, padding)
				; outputInt32 (c, doubleTag)
			)
			val dataSizePos = BinIO.getPosOut c
			val () = outputInt32 (c, 0)
		in
			{
				outstream = c
				, wholeSizePos = wholeSizePos
				, paddedSize = paddedSize
				, arraySizePos = arraySizePos
				, dataSizePos = dataSizePos
				, elements = ref 0
			}
		end
	fun writeDouble (da, f) = (
		outputReal (#outstream da, f)
		; #elements da := !(#elements da) + 1
	)
	fun concludeDoubleArray da =
		let
			val saved = BinIO.getPosOut (#outstream da)
		in
			BinIO.setPosOut (#outstream da, #wholeSizePos da)
			; outputInt32 (
				#outstream da
				, doubleArraySize (#paddedSize da, !(#elements da))
			); BinIO.setPosOut (#outstream da, #arraySizePos da)
			; outputInt32 (#outstream da, !(#elements da))
			; BinIO.setPosOut (#outstream da, #dataSizePos da)
			; outputInt32 (#outstream da, (!(#elements da) * 8))
			; BinIO.setPosOut (#outstream da, saved)
		end
	fun writeDoubleArray (c, name, array) =
		let
			val da = beginDoubleArray (c, name)
		in
			Vector.app (fn x => writeDouble (da, x)) array
			; concludeDoubleArray da
		end
	fun writeHeader (c, software, version) =
		let
			val date = Date.fromTimeUniv (Time.now ())
			val text = concat [
				"MATLAB 5.0 MAT-file, written by "
				, software, " ", version, ", "
				, Date.fmt "%Y-%m-$d %H:%M:%S UTC" date
			]
			val size = size text
			val header =
				CharVector.tabulate (
					124
					, fn i =>
						if i < size then String.sub (text, i)
						else #" "
				)
		in
			BinIO.output (c, Byte.stringToBytes header)
			; outputInt16 (c, 0x0100)
			; outputInt16 (c, 0x4d49)
		end
	fun openOut name =
		let
			val c = BinIO.openOut name
		in
			writeHeader (c, Program.name, Program.version)
			; c
		end
	val closeOut = BinIO.closeOut
end
