signature GZIP = sig
	exception Failure
	val openIn: string -> TextIO.instream
	val openOut: string -> TextIO.outstream
end

structure Gzip :> GZIP = struct
	structure Primitive :> sig
		eqtype gzFile
		val null: gzFile
		val gzopen: string * string -> gzFile
		val gzread: gzFile * CharArray.array * word * word -> int
		val gzwritea: gzFile * CharArray.array * word * word -> int
		val gzwritev: gzFile * CharVector.vector * word * word -> int
		val gzclose: gzFile -> int
	end = struct
		type gzFile = MLton.Pointer.t
		val null = MLton.Pointer.null
		val gzopen = _import "gzopen": string * string -> gzFile;
		val gzread = _import "gzreadoffset":
			gzFile * CharArray.array * word * word -> int;
		val gzwritea = _import "gzwriteoffset":
			gzFile * CharArray.array * word * word -> int;
		val gzwritev = _import "gzwriteoffset":
			gzFile * CharVector.vector * word * word -> int;
		val gzclose = _import "gzclose": gzFile -> int;
	end

	exception Failure
	fun readArraySlice (g, slice) =
		let
			val (array, offset, length) = CharArraySlice.base slice
			val wordoffset = Word.fromInt offset
			val wordlength = Word.fromInt length
		in
			Primitive.gzread (g, array, wordoffset, wordlength)
		end
	fun readerFromPrimitive (name, g) =
		let
			val closed = ref false
			fun error x = raise IO.Io {
				name = name
				, function = "readArr"
				, cause = x
			}
			fun readArr slice =
				if !closed then error IO.ClosedStream
				else let
					val r = readArraySlice (g, slice)
				in
					if r < 0 then error Failure
					else r
				end
			fun close () =
				if !closed then ()
				else (
					closed := true
					; ignore (Primitive.gzclose g)
				)
		in
			TextPrimIO.augmentReader (TextPrimIO.RD {
				name = name
				, chunkSize = 32 * 1024
				, readVec = NONE
				, readArr = SOME readArr
				, readVecNB = NONE
				, readArrNB = NONE
				, block = NONE
				, canInput = NONE
				, avail = fn () => NONE
				, getPos = NONE
				, setPos = NONE
				, endPos = NONE
				, verifyPos = NONE
				, close = close
				, ioDesc = NONE
			})
		end
	fun readerFromName name =
		let
			val path = name ^ "\000"
			val g = Primitive.gzopen (path, "r\000")
		in
			if g = Primitive.null then NONE
			else SOME (readerFromPrimitive (name, g))
		end
	fun openIn name =
		case readerFromName name of
			SOME x => TextIO.mkInstream (TextIO.StreamIO.mkInstream (x, ""))
			| NONE => raise IO.Io {
				name = name
				, function = "openIn"
				, cause = Failure
			}
	fun writeSlice (base, write) (g, slice) =
		let
			val (vector, offset, length) = base slice
			val wordoffset = Word.fromInt offset
			val wordlength = Word.fromInt length
		in
			write (g, vector, wordoffset, wordlength)
		end
	val writeVectorSlice = writeSlice (CharVectorSlice.base, Primitive.gzwritev)
	val writeArraySlice = writeSlice (CharArraySlice.base, Primitive.gzwritea)
	fun writerFromPrimitive (name, g) =
		let
			val closed = ref false
			fun error (function, x) = raise IO.Io {
				name = name
				, function = function
				, cause = x
			}
			fun close () =
				if !closed then ()
				else (
					closed := true
					; ignore (Primitive.gzclose g)
				)
			fun write (name, realWrite) slice =
				if !closed then error (name, IO.ClosedStream)
				else let
					val r = realWrite (g, slice)
				in
					if r <= 0 then error (name, Failure)
					else r
				end
			val writeVec = write ("writeVec", writeVectorSlice)
			val writeArr = write ("writeArr", writeArraySlice)
		in
			TextPrimIO.augmentWriter (TextPrimIO.WR {
				name = name
				, chunkSize = 32 * 1024
				, writeVec = SOME writeVec
				, writeArr = SOME writeArr
				, writeVecNB = NONE
				, writeArrNB = NONE
				, block = NONE
				, canOutput = NONE
				, getPos = NONE
				, setPos = NONE
				, endPos = NONE
				, verifyPos = NONE
				, close = close
				, ioDesc = NONE
			})
		end
	fun writerFromName name =
		let
			val path = name ^ "\000"
			val g = Primitive.gzopen (path, "w9")
		in
			if g = Primitive.null then NONE
			else SOME (writerFromPrimitive (name, g))
		end
	fun openOut name =
		case writerFromName name of
			SOME x => TextIO.mkOutstream (TextIO.StreamIO.mkOutstream (x, IO.NO_BUF))
			| NONE => raise IO.Io {
				name = name
				, function = "openOut"
				, cause = Failure
			}
end

