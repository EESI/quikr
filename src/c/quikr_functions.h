// our malloc checker
void check_malloc(void *ptr, char *error);

// count the sequences in a fasta file
unsigned long long count_sequences(const char *filename);

// normalize a matrix by dividing each element by the sum of it's column
void normalize_matrix(double *matrix, int height, int width);

// load a sensing matrix  
struct matrix *load_sensing_matrix(const char *filename, unsigned int target_kmer);

// get_rare_value 
void get_rare_value(double *count_matrix, unsigned long long width, double rare_percent, unsigned long long *ret_rare_value, unsigned long long  *ret_rare_width);
