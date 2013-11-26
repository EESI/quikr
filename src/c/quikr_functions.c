#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "kmer_utils.h"
#include "quikr.h"

void check_malloc(void *ptr, char *error) {
	if (ptr == NULL) {
		if(error != NULL)  {
			fprintf(stderr,"Error: %s\n", error); 
		} 
		else { fprintf(stderr, "Error: Could not allocate enough memory - %s\n", strerror(errno));
		} 
		exit(EXIT_FAILURE);
	}
}


void debug_arrays(double *count_matrix, struct matrix *sensing_matrix) {
	FILE *count_fh = fopen("count.mat", "w");
	FILE *sensing_fh = fopen("sensing.mat", "w");

	unsigned long long width = pow_four(sensing_matrix->kmer);
	unsigned long long i = 0;
	unsigned long long j = 0;

	for(i = 0; i < sensing_matrix->sequences; i++) {
		for(j = 0; j < width - 1; j++)
			fprintf(sensing_fh, "%lf\t", sensing_matrix->matrix[width*i + j]);
		fprintf(sensing_fh, "%lf\n", sensing_matrix->matrix[width*i + width-1]);
	}

	fclose(sensing_fh);

	for(j = 0; j < width - 1; j++)
		fprintf(count_fh, "%lf\t", count_matrix[j]);
	fprintf(count_fh, "%lf\n", count_matrix[width - 1]);

	fclose(count_fh);
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

unsigned long long count_sequences(const char *filename) { 
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  unsigned long long sequences = 0;

  FILE *fh = fopen(filename, "r");
  if(fh == NULL) {
    fprintf(stderr, "could not open \"%s\"\n", filename );
		return 0;
  }

  while ((read = getline(&line, &len, fh)) != -1) {
    if(line[0] == '>') 
      sequences++;
  }

  free(line);
  fclose(fh);

  return sequences;
}


struct matrix *load_sensing_matrix(const char *filename, unsigned int target_kmer) {

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
	check_malloc(line, NULL);

	// Check for quikr
	line = gzgets(fh, line, 1024);
	if(strcmp(line, "quikr\n") != 0) {
		fprintf(stderr, "This does not look like a quikr sensing matrix. Please check your path: %s\n", filename);
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
	if(kmer != target_kmer) {
		fprintf(stderr, "The sensing_matrix was trained with a different kmer than your requested kmer\n");
		exit(EXIT_FAILURE);
	}

	width = pow_four(kmer);

	// allocate a +1 size for the extra row
  matrix = malloc(sequences * (width) * sizeof(double));
	check_malloc(matrix, NULL);

	row = malloc((width) * sizeof(unsigned long long));
	check_malloc(row, NULL);
	
  headers = malloc(sequences * sizeof(char *));
	check_malloc(headers, NULL);

	for(i = 0; i < sequences; i++) {
		unsigned long long j = 0;
		// get header and add it to headers array 
		char *header = malloc(256 * sizeof(char));
		check_malloc(header, NULL);
		gzgets(fh, header, 256);
		if(header[0] != '>') {
			fprintf(stderr, "Error parsing sensing matrix, could not read header\n");
			exit(EXIT_FAILURE);
		}

		header[strlen(header) - 1] = '\0';
		headers[i] = header+1;

		row = memset(row, 0, (width) * sizeof(unsigned long long));

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

		}
		for(j = 0; j < width; j++) { 
			matrix[i*(width) + j] = ((double)row[j]);
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
