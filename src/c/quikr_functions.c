#include <stdio.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

int count_sequences(char *filename) { 
  char command[512];
  int sequences = 0;
  FILE *grep_output;

  sprintf(command, "grep -c ^\\> %s", filename);
  grep_output = popen(command, "r");
  if(grep_output == NULL) {
    fprintf(stderr, "Could not execute %s\n", command);
    exit(EXIT_FAILURE);
  }
  
  fscanf(grep_output, "%d", &sequences);

  pclose(grep_output);
  return sequences;
}


void normalize_matrix(double *matrix, int height, int width) {
  int x = 0;
  int y = 0;
  for(x = 0; x < height; x++) {

    double row_sum = 0;

    for(y = 0; y < (width); y++)
      row_sum = row_sum + matrix[width * x + y];
    for(y = 0; y < (width); y++) 
      matrix[width * x + y] = matrix[width * x + y] / row_sum; 
  }
}


double *load_count_matrix(char *filename, int width, int kmer) {

  double *count_matrix = malloc((width)*sizeof(double));
  char count_command[512];
  int x = 0;
  char *line = NULL;
  size_t len = 0;

  if(count_matrix == NULL) {
    fprintf(stderr, "could not allocate enough memory for the count matrix (%d x double) \n", width);
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


double *load_sensing_matrix(char *filename, int height, int width) {

  int x = 0;
  int y = 0;

  gzFile sensing_matrix_fh = NULL;

  double *sensing_matrix = malloc(height * width * sizeof(double));
  if(sensing_matrix == NULL) {
    fprintf(stderr, "Could not allocate memory for the sensing matrix\n");
  }

  sensing_matrix_fh = gzopen(filename, "r");
  if(sensing_matrix_fh == NULL) {
    fprintf(stderr, "could not open %s", filename);
    exit(EXIT_FAILURE);
  }

  // read our sensing matrix in
  for(x = 0; x < height; x++) {
    for (y = 1; y < width; y++) {
      char buffer[14];
      gzgets(sensing_matrix_fh, buffer, 14); 
      sensing_matrix[width*x + y] = atof(buffer);
    }
  }

  gzclose(sensing_matrix_fh);

  return sensing_matrix;
}

char **load_headers(char *filename, int sequences) {
  char command[512];
  char *line= NULL;
  int x = 0;
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

