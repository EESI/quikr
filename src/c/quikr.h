#define MATRIX_REVISION 0
#define pow_four(x) ( (unsigned long long)1 << (x * 2 ) )

struct matrix {
	unsigned long long sequences;
	unsigned int kmer;
	double *matrix;
	char **headers;
};
