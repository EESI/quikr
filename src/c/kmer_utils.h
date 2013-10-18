// Kmer functions
void convert_kmer_to_num(char *str, const unsigned long length);
unsigned long long * get_kmer_counts_from_file(const char *fn, const int kmer);
unsigned long num_to_index(const char *str, const int kmer, int *position, const long error_pos);

// Utility functions
char *strnstrip(const char *s, char *dest, int c, int len);
unsigned long long pow_four(unsigned long long x);

// Variables
const unsigned char alpha[256]; 
