#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "nnls.h"
#include "kmer_utils.h"
#include "quikr_functions.h"
#include "quikr.h"

#define sensing_matrix(i,j) (sensing_matrix[width*i + j])
#define USAGE "Usage:\n\tquikr [OPTION...] - Calculate estimated frequencies of bacteria in a sample.\n\nOptions:\n\n-i, --input\n\tthe sample's fasta file of NGS READS (fasta format)\n\n-s, --sensing-matrix\n\t location of the sensing matrix. (trained from quikr_train)\n\n-k, --kmer\n\tspecify what size of kmer to use. (default value is 6)\n\n-l, --lambda\n\tlambda value to use. (default value is 10000)\n\n-o, --output\n\tOTU_FRACTION_PRESENT a vector representing the percentage of database sequence's presence in sample. (csv output)\n\n-v, --verbose\n\tverbose mode.\n\n-V, --version\n\tprint version."

int main(int argc, char **argv) {

  int c;

  char *input_fasta_filename = NULL;
  char *sensing_matrix_filename = NULL;
  char *output_filename = NULL;

  unsigned long long x = 0;
  unsigned long long y = 0;
  unsigned long long z = 0;

  unsigned long long rare_value = 0;
  unsigned long long rare_width = 0; 

  double percentage = 1.0;
  double rare_percent = 1.0;

  unsigned long long width = 0;

  unsigned int kmer = 6;
  unsigned long long lambda = 10000;

  int verbose = 0;

  while (1) {
    static struct option long_options[] = {
      {"input", required_argument, 0, 'i'},
      {"kmer",  required_argument, 0, 'k'},
      {"lambda",  required_argument, 0, 'l'},
      {"output", required_argument, 0, 'o'},
      {"sensing-matrix", required_argument, 0, 's'},
      {"rare-percent", required_argument, 0, 'r'},
      {"verbose", no_argument, 0, 'v'},
      {"version", no_argument, 0, 'V'},
      {"help", no_argument, 0, 'h'},
      {"debug", no_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    int option_index = 0;

    c = getopt_long (argc, argv, "k:l:s:r:i:o:r:hdvV", long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {
      case 'k':
        kmer = atoi(optarg);
        break;
      case 'l':
        lambda = atoi(optarg);
        break;
      case 'r':
        rare_percent = atof(optarg);
        break;
      case 's':
        sensing_matrix_filename = optarg;
        break;
      case 'i':
        input_fasta_filename = optarg;
        break;
      case 'o':
        output_filename = optarg;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'V':
        printf("%s\n", VERSION);
        exit(EXIT_SUCCESS);
      case 'h':
        printf("%s\n", USAGE);
        exit(EXIT_SUCCESS);
      default:
        break;
    }
  }

  if(sensing_matrix_filename == NULL) {
    fprintf(stderr, "Error: sensing matrix filename (-s) must be specified\n\n");
    fprintf(stderr, "%s\n", USAGE);
    exit(EXIT_FAILURE);
  }
  if(output_filename == NULL) {
    fprintf(stderr, "Error: output filename (-o) must be specified\n\n");
    fprintf(stderr, "%s\n", USAGE);
    exit(EXIT_FAILURE);
  }
  if(input_fasta_filename == NULL) {
    fprintf(stderr, "Error: input fasta file (-i) must be specified\n\n");
    fprintf(stderr, "%s\n", USAGE);
    exit(EXIT_FAILURE);
  }

	if(rare_percent <= 0 || rare_percent > 1.0) {
    fprintf(stderr, "Error: rare percent must be between 0 and 1\n");
    exit(EXIT_FAILURE);
	}

  if(verbose) {
    printf("rare width:%ld\n", rare_width);
    printf("kmer: %u\n", kmer);
    printf("lambda: %llu\n", lambda);
    printf("fasta: %s\n", input_fasta_filename);
    printf("sensing matrix: %s\n", sensing_matrix_filename);
    printf("output: %s\n", output_filename);
  }

  if(access (sensing_matrix_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", sensing_matrix_filename);
    exit(EXIT_FAILURE);
  }
  if(access (input_fasta_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", input_fasta_filename);
    exit(EXIT_FAILURE);
  }

	if(kmer == 0) {
    fprintf(stderr, "Error: zero is not a valid kmer\n");
    exit(EXIT_FAILURE);
	}


  // 4 "ACGT" ^ Kmer gives us the size of output rows
  width = pow_four(kmer);

	// load counts matrix
	double *count_matrix = malloc(width * sizeof(double));
	check_malloc(count_matrix, NULL);

	{
		unsigned long long *integer_counts = get_kmer_counts_from_file(input_fasta_filename, kmer);

		for(x = 0; x < width; x++)
			count_matrix[x] = (double)integer_counts[x];

		free(integer_counts);
	}

	// load sensing matrix
	struct matrix *sensing_matrix = load_sensing_matrix(sensing_matrix_filename, kmer);

	// get our "rare" counts
	while(1) {
		rare_width = 0;
		for(x = 0; x < width; x++) {
			if(count_matrix[x] < rare_value) {
				rare_width++;
			}
		}
		percentage = (float)rare_width / (float)width;

		if(percentage >= rare_percent)
			break;

		rare_value++;
	}

	if(verbose)
		printf("there are %llu values less than %llu\n", rare_width, rare_value);

	// add a extra space for our zero's array
	rare_width++;

	// store our count matrix
	double *count_matrix_rare = calloc(rare_width, sizeof(double));
	check_malloc(count_matrix_rare, NULL);

	double *sensing_matrix_rare = malloc((rare_width) * sensing_matrix->sequences * sizeof(double));
	check_malloc(sensing_matrix_rare, NULL);

	// copy only kmers from our original counts that match our rareness percentage
	// in both our count matrix and our sensing matrix
	for(x = 0, y = 1;  x < width; x++) {
		if(count_matrix[x] < rare_value) {
			count_matrix_rare[y] = count_matrix[x];

			for(z = 0; z < sensing_matrix->sequences; z++) 
				sensing_matrix_rare[z*rare_width + y] = sensing_matrix->matrix[z*width + x];

			y++;
		}
	}

	// normalize our kmer counts and our sensing_matrix
	normalize_matrix(count_matrix_rare, 1, rare_width);
	normalize_matrix(sensing_matrix_rare, sensing_matrix->sequences, rare_width);

	// multiply our kmer counts and sensing matrix by lambda
	for(x = 1; x < rare_width; x++) 
		count_matrix_rare[x] *= lambda;

	for(x = 0; x < sensing_matrix->sequences; x++) {
		for(y = 1; y < rare_width; y++) {
			sensing_matrix_rare[rare_width*x + y] *= lambda;
		}
	}

	// count_matrix's first element should be zero
	count_matrix_rare[0] = 0;
	// stack one's on our first row of our sensing matrix
	for(x = 0; x < sensing_matrix->sequences; x++) {
		sensing_matrix_rare[x*rare_width] = 1.0;
	}

	// run NNLS
  double *solution = nnls(sensing_matrix_rare, count_matrix_rare, sensing_matrix->sequences, rare_width);

  // normalize our solution vector
  normalize_matrix(solution, 1, sensing_matrix->sequences);

  // output our matrix
  FILE *output_fh = fopen(output_filename, "w");
  if(output_fh == NULL) { 
    fprintf(stderr, "Could not open %s for writing\n", output_filename);
    exit(EXIT_FAILURE);
  }

  for(x = 0; x < sensing_matrix->sequences; x++)
      fprintf(output_fh, "%.10lf\n", solution[x]);

  fclose(output_fh);

	free(solution);

	free(count_matrix);
	free(sensing_matrix);

	free(count_matrix_rare);
	free(sensing_matrix_rare);

  return EXIT_SUCCESS;
}
