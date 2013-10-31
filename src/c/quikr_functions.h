// count the sequences in a fasta file
unsigned long long count_sequences(const char *filename);

// normalize a matrix by dividing each element by the sum of it's column
void normalize_matrix(double *matrix, int height, int width);

// load a sensing matrix  
struct matrix *load_sensing_matrix(const char *filename);

// add header and normalize count matrix
double *setup_count_matrix(char *filename, unsigned long long kmer, unsigned long long lambda, unsigned long long width);

// add header and normalize sensing matrix
struct matrix *setup_sensing_matrix(char *filename, unsigned long long kmer, unsigned long long lambda, unsigned long long width);
