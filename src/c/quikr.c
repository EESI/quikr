#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "nnls.h"
#include "quikr_functions.h"

#define sensing_matrix(i,j) (sensing_matrix[width*i + j])
#define USAGE "Usage:\n\tquikr [OPTION...] - Calculate estimated frequencies of bacteria in a sample.\n\nOptions:\n\n-i, --input\n\tthe sample's fasta file of NGS READS (fasta format)\n\n-f, --sensing-fasta\n\tlocation of the fasta file database used to create the sensing matrix (fasta format)\n\n-s, --sensing-matrix\n\t location of the sensing matrix. (trained from quikr_train)\n\n-k, --kmer\n\tspecify what size of kmer to use. (default value is 6)\n\n-l, --lambda\n\tlambda value to use. (default value is 10000)\n\n-o, --output\n\tOTU_FRACTION_PRESENT a vector representing the percentage of database sequence's presence in sample. (csv output)\n\n-v, --verbose\n\tverbose mode.\n\n-d, --debug\n\tdebug mode, read manpage for more details"

int main(int argc, char **argv) {


  int c;

  char *input_fasta_filename = NULL;
  char *sensing_matrix_filename = NULL;
  char *sensing_fasta_filename = NULL;
  char *output_filename = NULL;

  long x = 0;
  long y = 0;

  int verbose = 0;
  int debug = 0;

  int lambda = 10000;
  int kmer = 6;

  long width = 0;
  long sequences = 0;

  while (1) {
    static struct option long_options[] = {
      {"input", required_argument, 0, 'i'},
      {"kmer",  required_argument, 0, 'k'},
      {"lambda",  required_argument, 0, 'l'},
      {"output", required_argument, 0, 'o'},
      {"sensing-fasta",  required_argument, 0, 'f'},
      {"sensing-matrix", required_argument, 0, 's'},
      {"verbose", no_argument, 0, 'v'},
      {"debug", no_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    int option_index = 0;

    c = getopt_long (argc, argv, "k:l:f:s:i:o:hdv", long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {
      case 0:
      case 'k':
        kmer = atoi(optarg);
        break;
      case 'l':
        lambda = atoi(optarg);
        break;
      case 'f':
        sensing_fasta_filename = optarg;
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
      case 'd':
        debug = 1;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        printf("%s\n", USAGE);
        exit(EXIT_SUCCESS);
        break;
      default:
        break;
    }
  }


  if(sensing_matrix_filename == NULL) {
    fprintf(stderr, "Error: sensing matrix filename (-s) must be specified\n\n");
    fprintf(stderr, "%s\n", USAGE);
    exit(EXIT_FAILURE);
  }
  if(sensing_fasta_filename == NULL) {
    fprintf(stderr, "Error: sensing fasta filename (-f) must be specified\n\n");
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

  if(verbose) { 
    printf("kmer: %d\n", kmer);
    printf("lambda: %d\n", lambda);
    printf("fasta: %s\n", input_fasta_filename);
    printf("sensing database: %s\n", sensing_matrix_filename);
    printf("sensing database fasta: %s\n", sensing_fasta_filename);
    printf("output: %s\n", output_filename);
  }

  if(access (sensing_matrix_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", sensing_matrix_filename);
    exit(EXIT_FAILURE);
  }
  if(access (sensing_fasta_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", sensing_fasta_filename);
    exit(EXIT_FAILURE);
  }
  if(access (input_fasta_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", input_fasta_filename);
    exit(EXIT_FAILURE);
  }

  // 4 "ACGT" ^ Kmer gives us the size of output rows
  width = pow(4, kmer);
  width = width + 1;

  sequences = count_sequences(sensing_fasta_filename);
  if(sequences == 0) {
    fprintf(stderr, "Error: %s contains 0 fasta sequences\n", sensing_fasta_filename);
  }

  if(verbose) {
    printf("width: %ld\nsequences %ld\n", width, sequences);
  }

  double *sensing_matrix = load_sensing_matrix(sensing_matrix_filename, sequences, width);
  double *count_matrix = load_count_matrix(input_fasta_filename, width, kmer);

  // multiply our matrix by lambda
  for(x = 1; x < sequences; x++) {
    for(y= 0; y < width - 1; y++) {
      sensing_matrix(x, y) = sensing_matrix(x, y) * lambda;
    }
  }

  for(x = 0; x < sequences; x++) {
    sensing_matrix(x, 0) = 1.0;
  }
  // normalize our count_matrix
  normalize_matrix(count_matrix, 1, width);
  for(x = 0; x < width; x++) 
    count_matrix[x] = count_matrix[x] * lambda;
  
  // output our matricies if we are in verbose mode
  if(debug) { 
    FILE *sensing_matrix_fh = fopen( "sensing.matrix", "w");
    if(sensing_matrix_fh == NULL) {
      fprintf(stderr, "could not open sensing.matrix for writing.\n");
      exit(EXIT_FAILURE);
    }
    for(x = 0; x < sequences; x++) {
      for( y = 0; y < width; y++) {
        fprintf(sensing_matrix_fh, "%.10f\t", sensing_matrix(x, y));
      }
      fprintf(sensing_matrix_fh, "\n");
    }
    fclose(sensing_matrix_fh);

    FILE *count_matrix_fh = fopen("count.matrix", "w");
    if(count_matrix_fh == NULL) {
      fprintf(stderr, "could not open sensing.matrix for writing.\n");
      exit(EXIT_FAILURE);
    }
    for(x = 0; x < width; x++) {
      fprintf(count_matrix_fh, "%.10f\n", count_matrix[x]);
    }
    fclose(count_matrix_fh);
  }

  double *solution = nnls(sensing_matrix, count_matrix, sequences, width);

  // normalize our solution vector
  normalize_matrix(solution, 1, sequences);

  // output our matrix
  FILE *output_fh = fopen(output_filename, "w");
  if(output_fh == NULL) { 
    fprintf(stderr, "Could not open %s for writing\n", output_filename);
    exit(EXIT_FAILURE);
  }
  for(x = 0; x < sequences; x++) {
      fprintf(output_fh, "%.10lf\n", solution[x]);
  }
  fclose(output_fh);

  return EXIT_SUCCESS;
}
