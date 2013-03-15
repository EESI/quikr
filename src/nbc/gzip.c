#include <zlib.h>

int
gzreadoffset(gzFile file, voidp buf, unsigned offset, unsigned len)
{
	return gzread(file, buf + offset, len);
}

int
gzwriteoffset(gzFile file, voidpc buf, unsigned offset, unsigned len)
{
	return gzwrite(file, buf + offset, len);
}
