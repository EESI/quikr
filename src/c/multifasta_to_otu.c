#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "nnls.h"
#include "quikr_functions.h"

#define sensing_matrix(i,j) (sensing_matrix[width*i + j])
#define solutions(i,j) (solutions[sequences*i+ j])
#define USAGE "Usage:\n\tmultifasta_to_otu [OPTION...] - create a QIIME OTU table based on Quikr results. \n\nOptions:\n\n-i, --input-directory\n\tthe directory containing the samples' fasta files of reads (note each file should correspond to a separate sample)\n\n-f, --sensing-fasta\n\tlocation of the fasta file database used to create the sensing matrix (fasta format)\n\n-s, --sensing-matrix\n\t location of the sensing matrix. (sensing from quikr_train)\n\n-k, --kmer\n\tspecify what size of kmer to use. (default value is 6)\n\n-l, --lambda\n\tlambda value to use. (default value is 10000)\n\n-j, --jobs\n\t specifies how many jobs to run at once. (default value is the number of CPUs)\n\n-o, --output\n\tthe OTU table, with NUM_READS_PRESENT for each sample which is compatible with QIIME's convert_biom.py (or a sequence table if not OTU's)\n\n-v, --verbose\n\tverbose mode.\n\n-V, --version\n\tprint version."

int main(int argc, char **argv) {

  int c;

  char *input_fasta_directory = NULL;
  char *sensing_matrix_filename = NULL;
  char *sensing_fasta_filename = NULL;
  char *output_filename = NULL;

  double *sensing_matrix;

  long width = 0;
  long sequences = 0;

  int kmer = 6;
  int lambda = 10000;

  long x = 0;
  long y = 0;

  int jobs = 1;
  #ifdef Linux
  jobs = get_nprocs();
  #endif
  #ifdef Darwin
  jobs = sysconf (_SC_NPROCESSORS_ONLN);
  #endif

  int verbose = 0;

  DIR *input_directory_dh;
  struct dirent *entry;

  while (1) {
    static struct option long_options[] = {
      {"input-directory", required_argument, 0, 'i'},
      {"kmer",  required_argument, 0, 'k'},
      {"lambda",  required_argument, 0, 'l'},
      {"jobs",  required_argument, 0, 'j'},
      {"output", required_argument, 0, 'o'},
      {"sensing-fasta",  required_argument, 0, 'f'},
      {"sensing-matrix", required_argument, 0, 's'},
      {"verbose", no_argument, 0, 'v'},
      {"version", no_argument, 0, 'V'},
      {0, 0, 0, 0}
    };
    int option_index = 0;

    c = getopt_long (argc, argv, "k:l:f:s:i:o:j:hvV", long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {
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
      case 'j':
        jobs = atoi(optarg);
        break;
      case 'i':
        input_fasta_directory = optarg;
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
        break;
      case 'h':
        puts(USAGE);
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
  if(input_fasta_directory == NULL) {
    fprintf(stderr, "Error: input fasta directory (-i) must be specified\n\n");
    fprintf(stderr, "%s\n", USAGE);
    exit(EXIT_FAILURE);
  }

  if(verbose) { 
    printf("kmer: %d\n", kmer);
    printf("lambda: %d\n", lambda);
    printf("input directory: %s\n", input_fasta_directory);
    printf("sensing database: %s\n", sensing_matrix_filename);
    printf("sensing database fasta: %s\n", sensing_fasta_filename);
    printf("output: %s\n", output_filename);
    printf("number of jobs to run at once: %d\n", jobs); 
  }


  if(access (sensing_matrix_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", sensing_matrix_filename);
    exit(EXIT_FAILURE);
  }

  if(access (sensing_fasta_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", sensing_fasta_filename);
    exit(EXIT_FAILURE);
  }

  input_directory_dh = opendir(input_fasta_directory);
  if(input_directory_dh == NULL) {
    fprintf(stderr, "could not open %s\n", input_fasta_directory);
    exit(EXIT_FAILURE);
  } 

  // do a directory count
  long dir_count = -2; // -2 for ../ and ./
  while(entry = readdir(input_directory_dh)) 
    dir_count++; 
  rewinddir(input_directory_dh);
  if(dir_count == 0) {
    fprintf(stderr, "%s is empty\n", input_fasta_directory);
    exit(EXIT_FAILURE);
  }

  // 4 "ACGT" ^ Kmer gives us the size of output rows
  width = pow(4, kmer) + 1;
  sequences = count_sequences(sensing_fasta_filename);
  if(sequences == 0) {
    fprintf(stderr, "Error: %s contains 0 fasta sequences\n", sensing_fasta_filename);
  }

  if(verbose) {
    printf("directory count: %ld\n", dir_count);
    printf("width: %ld\nsequences %ld\n", width, sequences);
  }

  sensing_matrix = load_sensing_matrix(sensing_matrix_filename, sequences, width);

  // multiply our matrix by lambda
  for(x = 0; x < sequences; x++) {
    for(y= 0; y < width; y++) {
      sensing_matrix(x, y) = sensing_matrix(x, y) * lambda;
    }
  }

  // set the first row to be all 1's
  for(x = 0; x < sequences; x++) {
    sensing_matrix(x, 0) = 1.0;
  }

  double *solutions = malloc(dir_count * sequences * sizeof(double));
  if(solutions == NULL) {
    fprintf(stderr, "Could not allocate enough memory for solutions vector\n");
    exit(EXIT_FAILURE);
  }

  char **filenames = malloc(dir_count * sizeof(char *));
  if(filenames == NULL) {
    fprintf(stderr, "Could not allocate enough memory\n");
    exit(EXIT_FAILURE);
  }

  long *file_sequence_count = malloc(dir_count * sizeof(long));
  if(file_sequence_count == NULL) {
    fprintf(stderr, "Could not allocate enough memory\n");
    exit(EXIT_FAILURE);
  }

  omp_set_num_threads(jobs);
  long done = 0;
  printf("Beginning to process samples\n"); 
  #pragma omp parallel for shared(solutions, sequences, width, done)
  for(long i = 0; i < dir_count; i++ ) {

    long z = 0;
    struct dirent *directory_entry;
    char *filename = malloc(256 * sizeof(char));
    char *base_filename = malloc(256 * sizeof(char));
    if(filename == NULL || base_filename == NULL) {
      fprintf(stderr, "Could not allocate enough memory\n");
      exit(EXIT_FAILURE);
    }

    #pragma omp critical
    {
      directory_entry = readdir(input_directory_dh);
      strcpy(base_filename, directory_entry->d_name);
    }

    if(strcmp(base_filename, "..") == 0 || strcmp(base_filename, ".") == 0) {
      i--;
      continue;
    }

    // get our base filenames
    filenames[i] = base_filename; 

    // get our real filename
    sprintf(filename, "%s/%s", input_fasta_directory, directory_entry->d_name);

    // get individual sequence count
    file_sequence_count[i] =  count_sequences(filename);

    // count the kmer amounts 
    double *count_matrix = load_count_matrix(filename, width, kmer);

    // normalize our kmer counts
    normalize_matrix(count_matrix, 1, width);

    // multiply our kmers frequency by lambda
    for(z = 0; z < width; z++) 
      count_matrix[z] = count_matrix[z] * lambda;

    double *sensing_matrix_copy = malloc(sizeof(double) * sequences * width);
    if(sensing_matrix_copy == NULL) {
      fprintf(stderr, "Could not allocate enough memory\n");
      exit(EXIT_FAILURE);
    }

    memcpy(sensing_matrix_copy, sensing_matrix, sequences * width * sizeof(double));


    // run nnls
    double *solution = nnls(sensing_matrix_copy, count_matrix, sequences, width);

    // normalize our solution
    normalize_matrix(solution, 1, sequences);

    // add the current solution to the solutions array
    for(z = 0; z < sequences; z++ )  {
      solutions(i, z) = solution[z]; 
    }

    done++;
    printf("%ld/%ld samples processed\n", done, dir_count); 
    free(solution);
    free(count_matrix);
    free(filename);
    free(sensing_matrix_copy);
  }

  char **headers = load_headers(sensing_fasta_filename, sequences);

  // output our matrix
  FILE *output_fh = fopen(output_filename, "w");
  if(output_fh == NULL) { 
    fprintf(stderr, "Could not open %s for writing\n", output_filename);
    exit(EXIT_FAILURE);
  }

  fprintf(output_fh, "# QIIME vQuikr OTU table\n");
  fprintf(output_fh, "#OTU_ID\t");

  // print our filename headers
  for(x = 0; x < dir_count - 1; x++) {
    fprintf(output_fh, "%s\t", filenames[x]);
  }
  fprintf(output_fh, "%s\n", filenames[dir_count - 1]);

  // get our actual values 
  for(y = 0; y < sequences; y++) {
    for(x = 0; x < dir_count; x++) {
      solutions(x, y) = round(solutions(x, y) * file_sequence_count[x]);
    }
  }

  for(y = 0; y < sequences; y++) {

    double column_sum = 0.;
    for(x = 0; x < dir_count; x++) {
      column_sum += solutions(x, y);
    }

    // if our column is zero, don't bother printing the row
    if(column_sum != 0) {
      fprintf(output_fh, "%s\t", headers[y]);

      for(x = 0; x < dir_count - 1; x++) {
        fprintf(output_fh, "%d\t", (int)solutions(x, y));
      }
      fprintf(output_fh, "%d\n", (int)solutions[sequences*(dir_count - 1) + y]);
    }
  }
  fclose(output_fh);

  return EXIT_SUCCESS;
}
