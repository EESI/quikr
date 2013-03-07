# Helper Functions

# Check if file is gzip compressed
def isCompressed(filename):
  f = open(filename, "rb")

  # The first two bytes of a tar file are always '1f 8b'
  if f.read(2) == '\x1f\x8b':
    f.close()
    return True
  else:
    f.close()
    return False
 

