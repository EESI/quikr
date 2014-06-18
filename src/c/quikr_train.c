#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "kmer_utils.h"
#include "quikr_functions.h"

#define USAGE "Usage:\n\tquikr_train [OPTION...] - train a database for use with quikr.\n\nOptions:\n\n-i, --input\n\tthe database of sequences to create the sensing matrix (fasta format)\n\n-k, --kmer\n\tspecify what size of kmer to use. (default value is 6)\n\n-o, --output\n\tthe sensing matrix. (a gzip'd text file)\n\n-v, --verbose\n\tverbose mode.\n\n-V, --version\n\tprint version."

int main(int argc, char **argv) {

	// getline variables
  char *line = NULL;
  size_t len = 0;
  ssize_t read;


  int c;

	// k-mer is 6 by default
  int kmer = 6;

	// revision number
	int revision = 0;
	// iterators
  long long i = 0;
  unsigned long long j = 0;
	unsigned long long position = 0;

  int verbose = 0;
  int force_name = 0;

  char *fasta_filename = NULL;
  char *output_file = NULL;

  gzFile output = NULL;
  FILE *input = NULL;

  while (1) {
    static struct option long_options[] = {
      {"verbose", no_argument, 0, 'v'},
      {"force_name", no_argument, 0, 'f'},
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {"input", required_argument, 0, 'i'},
      {"kmer",  required_argument, 0, 'k'},
      {"output", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    int option_index = 0;

    c = getopt_long (argc, argv, "i:o:k:fhvV", long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {
      case 'i':
        fasta_filename = optarg;
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
      case 'f':
        force_name = 1;
        break;
      case 'V':
        printf("%s\n", VERSION);
        exit(EXIT_SUCCESS);
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

  if(fasta_filename == NULL) {
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
    printf("fasta file: %s\n", fasta_filename);
    printf("output file: %s\n", output_file);
  }

  if(access (fasta_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", fasta_filename);
    exit(EXIT_FAILURE);
  }

	if(kmer == 0) {
    fprintf(stderr, "Error: zero is not a valid kmer\n");
    exit(EXIT_FAILURE);
	}

  if(strcmp(&output_file[strlen(output_file) - 3], ".gz") != 0 && !force_name) {
    char *temp = malloc(strlen(output_file) + 4);
    if(temp == NULL) {
      fprintf(stderr, "Could not allocate enough memory\n"); 
      exit(EXIT_FAILURE);
    }
    sprintf(temp, "%s.gz", output_file);
    output_file = temp;
    printf("appending a .gz to our output file: %s\n", output_file);
  }

	// 4 ^ Kmer gives us the width, or the number of permutations of ACTG with
	// kmer length
  unsigned long width = pow(4, kmer);
  unsigned long long sequences = count_sequences(fasta_filename);
  if(sequences == 0) {
    fprintf(stderr, "Error: %s contains 0 fasta sequences\n", fasta_filename);
		exit(EXIT_FAILURE);
  }

  if(verbose) {
    printf("sequences: %llu\nwidth: %ld\n", sequences, width);
    printf("Writing our sensing matrix to %s\n", output_file);
	}

  input = fopen(fasta_filename, "r" );
  if(input == NULL) {
    fprintf(stderr, "Error opening %s - %s\n", fasta_filename, strerror(errno));
    exit(EXIT_FAILURE);
  }

  // open our output file
  output = gzopen(output_file, "w");
  if(output == NULL) {
    fprintf(stderr, "Error: could not open output file, error code: %s\n", strerror(errno));
    exit(EXIT_FAILURE);
  }

	// create our header 
	gzprintf(output, "quikr\n");
	gzprintf(output, "%ld\n", revision);
	gzprintf(output, "%ld\n", sequences);
	gzprintf(output, "%d\n", kmer);

	// malloc our return array
  unsigned long long * counts = malloc((width + 1) * sizeof(unsigned long long));
  if(counts == NULL)  {
		fprintf(stderr, strerror(errno));
    exit(EXIT_FAILURE);
	}

	char *str = malloc(4096);
	if(str == NULL) { 
		fprintf(stderr, strerror(errno));
		exit(EXIT_FAILURE);
	}

	unsigned long long str_size = 4096;

	// seek the first character, and skip over it
	fseek(input, 1, SEEK_CUR);
	while ((read = getseq(&line, &len, input)) != -1) {

		// find first whitespace
		for(i = 0; i < read; i ++) {
			if(line[i] == ' ' || line[i] == '\t' || line[i] == '\n')
				break;
		}
		
		// write our header
    gzprintf(output, ">%.*s\n", i, line);

		// find our first \n, this should be the end of the header
		char *start = strchr(line, '\n');	
		if(start == NULL) 
			continue;

		size_t start_len = strlen(start);


		// if our current str buffer isn't big enough, realloc
		if(start_len + 1 > str_size + 1) { 
			str = realloc(str, start_len + 1);
			if(str == NULL) { 
				exit(EXIT_FAILURE);
				fprintf(stderr, strerror(errno));
			}
		}

		// strip out all other newlines to handle multiline sequences
		str = strnstrip(start, str, '\n',start_len);
		size_t seq_length = strlen(str);

		// relace A, C, G and T with 0, 1, 2, 3 respectively
		// everything else is 5 
		for(j = 0; j < seq_length; j++) {
			str[j] = alpha[(int)str[j]];
		}
		
		// set counts to zero
		memset(counts, 0, width * sizeof(unsigned long long));

		// loop through our string to process each k-mer
		for(position = 0; position < (seq_length - kmer + 1); position++) {
			unsigned long mer = 0;
			unsigned long multiply = 1;

			// for each char in the k-mer check if it is an error char
			for(i = position + kmer - 1; i >= (signed)position; i--){
				if(str[i] >> 2) { 
					mer = width;
					position = i;
					goto next;
				}

				// multiply this char in the mer by the multiply
				// and bitshift the multiply for the next round
				mer += str[i] * multiply;
				multiply = multiply << 2;
			}
			// use this point to get mer of our loop
			next:
			// bump up the mer value in the counts array
			counts[mer]++;
		}

		for(j = 0; j < width; j++) {
			gzprintf(output, "%lld\n", counts[j]);
		}

	} 

	free(counts);
	free(line);

  gzclose(output);
  fclose(input);

  return EXIT_SUCCESS;
}
