#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "quikr_functions.h"

#define USAGE "Usage:\n\tquikr_train [OPTION...] - to train a database for use with quikr.\n\nOptions:\n\n-i, --input\n\tthe database of sequences to create the sensing matrix (fasta format)\n\n-k, --kmer\n\tspecify what size of kmer to use. (default value is 6)\n\n-o, --output\n\tthe sensing matrix. (a gzip'd text file)\n\n-v, --verbose\n\tverbose mode."

int main(int argc, char **argv) {

  char probabilities_command[512];
  char kmers_file[256];
  char *line = NULL;
  char *val;
  size_t len = 0;


  int c;
  int kmer = 6;

  char *fasta_file = NULL;
  char *output_file = NULL;

  int x = 0;
  int y = 0;

  int verbose = 0;

  gzFile output = NULL;

  while (1) {
    static struct option long_options[] = {
      {"verbose", no_argument, 0, 'v'},
      {"input", required_argument, 0, 'i'},
      {"kmer",  required_argument, 0, 'k'},
      {"output", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    int option_index = 0;

    c = getopt_long (argc, argv, "i:o:k:hv", long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {
      case 'i':
        fasta_file = optarg;
        break;
      case 'k':
        kmer = atoi(optarg);
        break;
      case 'o':
        output_file = optarg;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        printf("%s\n", USAGE);
        exit(EXIT_SUCCESS);
        break;
      case '?':
        /* getopt_long already printed an error message. */
        break;
      default:
        exit(EXIT_FAILURE);
    }
  }


  if(fasta_file == NULL) {
    fprintf(stderr, "Error: input fasta file (-i) must be specified\n\n");
    fprintf(stderr, "%s\n", USAGE);
    exit(EXIT_FAILURE);
  }

  if(output_file == NULL) {
    fprintf(stderr, "Error: output matrix file (-o) must be specified\n\n");
    fprintf(stderr, "%s\n", USAGE);
    exit(EXIT_FAILURE);
  }

  if(verbose) {
    printf("kmer size: %d\n", kmer);
    printf("fasta file: %s\n", fasta_file);
    printf("output file: %s\n", output_file);
  }

  if(access (fasta_file, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", fasta_file);
    exit(EXIT_FAILURE);
  }

  if(strcmp(&output_file[strlen(output_file) - 3], ".gz") != 0) {
    char *temp = malloc(strlen(output_file) + 4);
    if(temp == NULL) {
      fprintf(stderr, "Could not allocate enough memory\n"); 
      exit(EXIT_FAILURE);
    }
    sprintf(temp, "%s.gz", output_file);
    output_file = temp;
    printf("appending a .gz to our output file: %s\n", output_file);
  }

  // 4 ^ Kmer gives us the width, or the number of permutations of ACTG with kmer length
  int width = pow(4, kmer);
  int sequences = count_sequences(fasta_file);
  if(sequences == 0) {
    fprintf(stderr, "Error: %s contains 0 fasta sequences\n", fasta_file);
  }

  if(verbose)
    printf("sequences: %d\nwidth: %d\n", sequences, width);

  // Allocate our matrix with the appropriate size, just one row
  double *trained_matrix = malloc(width*sizeof(double));
  if(trained_matrix == NULL) {
    fprintf(stderr, "Could not allocate enough memory\n");
    exit(EXIT_FAILURE);
  }

  // call the probabilities-by-read command
  sprintf(probabilities_command, "generate_kmers %d | probabilities-by-read %d %s /dev/stdin", kmer, kmer, fasta_file);
  FILE *probabilities_output = popen(probabilities_command, "r");
  if(probabilities_output == NULL) {
    fprintf(stderr, "Error could not execute: %s\n", probabilities_command);
    exit(EXIT_FAILURE);
  }

  // open our output file
  output = gzopen(output_file, "w6");
  if(output == NULL) {
    fprintf(stderr, "Error: could not open output file, error code: %d", errno);
    exit(EXIT_FAILURE);
  }

  if(verbose)
    printf("Writing our sensing matrix to %s\n", output_file);

  // read normalize and write our matrix  in one go
  for(x = 0; x < sequences; x++) {

    getline(&line, &len, probabilities_output);

    // Read our first element in outside of the loop
    val = strtok(line,"\t\n\r");
    trained_matrix[0] = atof(val);
    // iterate through and load the array
    for (y = 1; y < width; y++) {
      val = strtok (NULL, "\t\n\r");
      trained_matrix[y] = atof(val);
    }

    double row_sum = 0;

    for( y = 0; y < width; y++) {
      row_sum = row_sum + trained_matrix[y];
    }

    for( y = 0; y < width; y++) {
      trained_matrix[y] = trained_matrix[y] / row_sum;
    }

    for( y = 0; y < width; y++) {
      gzprintf(output, "%.10f\t", trained_matrix[y]);
    }
    gzprintf(output, "\n");
  }

  free(trained_matrix);
  gzclose(output);
  pclose(probabilities_output);

  return 0;
}
