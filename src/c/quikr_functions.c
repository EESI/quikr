#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "quikr.h"

unsigned long long count_sequences(const char *filename) { 

  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  unsigned long long sequences = 0;

  FILE *fh = fopen(filename, "r");
  if(fh == NULL) {
    fprintf(stderr, "could not open\"%s\"", filename );
    exit(EXIT_FAILURE);
  }

  while ((read = getline(&line, &len, fh)) != -1) {
    if(line[0] == '>') 
      sequences++;
  }

  free(line);
  fclose(fh);

  return sequences;
}


void normalize_matrix(double *matrix, unsigned long long height, unsigned long long width) {
  unsigned long long x = 0;
  unsigned long long y = 0;

  for(x = 0; x < height; x++) {

    double row_sum = 0;

    for(y = 0; y < (width); y++)
      row_sum = row_sum + matrix[width * x + y];
    for(y = 0; y < (width); y++) 
      matrix[width * x + y] = matrix[width * x + y] / row_sum; 
  }
}


double *load_count_matrix(const char *filename, const unsigned long long width, const unsigned int kmer) {

  double *count_matrix = malloc(width*sizeof(double));
  char count_command[1024];
  unsigned long long x = 0;
  char *line = NULL;
  size_t len = 0;

  if(count_matrix == NULL) {
    fprintf(stderr, "could not allocate enough memory for the count matrix\n");
    exit(EXIT_FAILURE);
  }

  // create out count matrix
  sprintf(count_command, "count-kmers -r %d -1 -u %s", kmer, filename);
  FILE *count_output = popen(count_command, "r");
  if(count_output == NULL) {
    fprintf(stderr, "could not execute \"%s\"", count_command);
    exit(EXIT_FAILURE);
  }

  // set first element to zero.
  count_matrix[0] = 0;

  // get our first line
  getline(&line, &len, count_output);
  count_matrix[1] = atoi(line);

  // iterate over the rest of the lines
  for(x = 2; x < width; x++) { 
    getline(&line, &len, count_output);
    count_matrix[x] = atoi(line);
  } 

  free(line);
  pclose(count_output);

  return count_matrix;
}


struct matrix *load_sensing_matrix(const char *filename) {

	char *line = NULL;
	char **headers = NULL;

	double *matrix = NULL;

	int kmer = 0;
	
	unsigned long long i = 0;	
	unsigned long long *row = NULL;
	unsigned long long sequences = 0;
	unsigned long long width = 0;

 	struct matrix *ret = NULL;

  gzFile fh = NULL;

  fh = gzopen(filename, "r");
  if(fh == NULL) {
    fprintf(stderr, "could not open %s", filename);
    exit(EXIT_FAILURE);
  }

	line = malloc(1024 * sizeof(char));
	if(line == NULL) 
		exit(EXIT_FAILURE);

	// Check for quikr
	line = gzgets(fh, line, 1024);
	if(strcmp(line, "quikr\n") != 0) {
		fprintf(stderr, "This does not look like a quikr sensing matrix. Please check your path\n");
		exit(EXIT_FAILURE);
	}

	// check version
	line = gzgets(fh, line, 1024);
	if(atoi(line) != MATRIX_REVISION) {
		fprintf(stderr, "Sensing Matrix uses an unsupported version, please retrain your matrix\n");
		exit(EXIT_FAILURE);
	}

	// get number of sequences
	line = gzgets(fh, line, 1024);
	sequences = strtoull(line, NULL, 10);
	if(sequences == 0) {
		fprintf(stderr, "Error parsing sensing matrix, sequence count is zero\n");
		exit(EXIT_FAILURE);
	}

	// get kmer
	gzgets(fh, line, 1024);
	kmer = atoi(line);
	if(kmer == 0) {
		fprintf(stderr, "Error parsing sensing matrix, kmer is zero\n");
		exit(EXIT_FAILURE);
	}

	width = pow_four(kmer);

	// allocate a +1 size for the extra row
  matrix = malloc(sequences * (width + 1) * sizeof(double));
  if(matrix == NULL) {
    fprintf(stderr, "Could not allocate memory for the sensing matrix\n");
  }

	row = malloc(width * sizeof(unsigned long long));
	if(row == NULL) {
		fprintf(stderr, "Could not allocate memory for parsing row\n");
	}
	
  headers = malloc(sequences * sizeof(char *));
  if(headers == NULL) {
    fprintf(stderr, "could not allocate enough memory for header pointers\n");
    exit(EXIT_FAILURE);
  }

	for(i = 0; i < sequences; i++) {
		unsigned long long sum = 0;
		unsigned long long j = 0;
		// get header and add it to headers array 
		char *header = malloc(256 * sizeof(char));
		gzgets(fh, header, 256);
		if(header[0] != '>') {
			fprintf(stderr, "Error parsing sensing matrix, could not read header\n");
			exit(EXIT_FAILURE);
		}

		headers[i] = header+1;

		row = memset(row, 0, (width + 1) * sizeof(unsigned long long));

		for(j = 0; j < width; j++) {
			line = gzgets(fh, line, 32);
			if(line == NULL || line[0] == '>') {
				fprintf(stderr, "Error parsing sensing matrix, line does not look like a value\n");
				exit(EXIT_FAILURE);
			}

			row[j] = strtoull(line, NULL, 10);
			if(errno) {
				printf("could not parse '%s'\n into a number", line);
				exit(EXIT_FAILURE);
			}

			sum += row[j];
		}
		for(j = 1; j < width+1; j++) { 
			matrix[i*(width+1) + j] = ((double)row[j-1]) / sum;
		}
	}

	// load the matrix of counts
  gzclose(fh);

	free(line);
	free(row);

	ret = malloc(sizeof(struct matrix));
	(*ret).kmer = kmer;
	(*ret).sequences = sequences;
	(*ret).matrix = matrix;
	(*ret).headers = headers;

  return ret;
}


char **load_headers(char *filename, long sequences) {
  char command[512];
  char *line= NULL;
  long x = 0;
  FILE *grep_output;
  size_t len = 0;

  sprintf(command, "grep ^\\> %s", filename);
  grep_output = popen(command, "r");
  if(grep_output == NULL) {
    fprintf(stderr, "Could not execute %s\n", command);
    exit(EXIT_FAILURE);
  }
  
  char **headers = malloc(sequences * sizeof(char *));
  if(headers == NULL) {
    fprintf(stderr, "could not allocated enough memory\n");
    exit(EXIT_FAILURE);
  }

  for(x = 0; x < sequences; x++) {

    char *header = malloc(256 * sizeof(char));
    if(header == NULL) {
      fprintf(stderr, "could not allocated enough memory\n");
      exit(EXIT_FAILURE);
    }
    getline(&line, &len, grep_output);
    sscanf(line + 1, "%s", header);
    headers[x] = header;
  }

  free(line);
  pclose(grep_output);

  return headers;
}

