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

#include "kmer_utils.h"
#include "nnls.h"
#include "quikr.h"
#include "quikr_functions.h"

#ifdef Linux
#include <sys/sysinfo.h>
#endif

#define sensing_matrix(i,j) (sensing_matrix[width*i + j])
#define solutions(i,j) (solutions[sequences*i+ j])
#define USAGE "Usage:\n\tmultifasta_to_otu [OPTION...] - create a QIIME OTU table based on Quikr results. \n\nOptions:\n\n-i, --input-directory\n\tthe directory containing the samples' fasta files of reads (note each file should correspond to a separate sample)\n\n-f, --sensing-fasta\n\tlocation of the fasta file database used to create the sensing matrix (fasta format)\n\n-s, --sensing-matrix\n\t location of the sensing matrix. (sensing from quikr_train)\n\n-k, --kmer\n\tspecify what size of kmer to use. (default value is 6)\n\n-l, --lambda\n\tlambda value to use. (default value is 10000)\n\n-j, --jobs\n\t specifies how many jobs to run at once. (default value is the number of CPUs)\n\n-o, --output\n\tthe OTU table, with NUM_READS_PRESENT for each sample which is compatible with QIIME's convert_biom.py (or a sequence table if not OTU's)\n\n-v, --verbose\n\tverbose mode.\n\n-V, --version\n\tprint version."

int main(int argc, char **argv) {

  int c;

  char *input_fasta_directory = NULL;
  char *sensing_matrix_filename = NULL;
  char *output_filename = NULL;

  unsigned long long x = 0;
  unsigned long long y = 0;

 	long long i = 0;

  unsigned long long width = 0;

  unsigned int kmer = 6;
  unsigned long long lambda = 10000;


  unsigned int jobs = 1;

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
      {"sensing-matrix", required_argument, 0, 's'},
      {"verbose", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {0, 0, 0, 0}
    };
    int option_index = 0;

    c = getopt_long (argc, argv, "k:l:s:i:o:j:hvV", long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {
      case 'k':
        kmer = atoi(optarg);
        break;
      case 'l':
        lambda = atoi(optarg);
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
    printf("kmer: %u\n", kmer);
    printf("lambda: %llu\n", lambda);
    printf("input directory: %s\n", input_fasta_directory);
    printf("sensing database: %s\n", sensing_matrix_filename);
    printf("output: %s\n", output_filename);
    printf("number of jobs to run at once: %d\n", jobs); 
  }


  if(access (sensing_matrix_filename, F_OK) == -1) {
    fprintf(stderr, "Error: could not find %s\n", sensing_matrix_filename);
    exit(EXIT_FAILURE);
  }

  input_directory_dh = opendir(input_fasta_directory);
  if(input_directory_dh == NULL) {
    fprintf(stderr, "could not open %s\n", input_fasta_directory);
    exit(EXIT_FAILURE);
  } 

	if(kmer == 0) {
    fprintf(stderr, "Error: zero is not a valid kmer\n");
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

  struct matrix *sensing_matrix = load_sensing_matrix(sensing_matrix_filename);

  if(verbose) {
    printf("directory count: %ld\n", dir_count);
    printf("width: %llu\nsequences %llu\n", width, sensing_matrix->sequences);
  }

  // multiply our matrix by lambda
  for(x = 0; x < sensing_matrix->sequences; x++) {
    for(y= 0; y < width; y++) {
      sensing_matrix->matrix[x*width + y] *= lambda;
    }
  }

  // set the first row to be all 1's
  for(x = 0; x < sensing_matrix->sequences; x++) {
    sensing_matrix->matrix[x*width] = 1.0;
  }

  double *solutions = malloc(dir_count * sensing_matrix->sequences * sizeof(double));
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
  #pragma omp parallel for shared(solutions, sensing_matrix, width, done)
  for(long i = 0; i < dir_count; i++ ) {

    unsigned long long z = 0;
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

    // count the kmer amounts, and convert it to a double array
		unsigned long long *integer_counts = get_kmer_counts_from_file(filename, kmer);
		double *count_matrix = malloc(sizeof(double) * width);
		if(count_matrix == NULL) {
			fprintf(stderr, "Could not allocate memory:\n");
			exit(EXIT_FAILURE);
		}

		count_matrix[0] = 0; 

		for(x = 0; x < width - 1 ; x++)
			count_matrix[x+1] = (double)integer_counts[x];

		free(integer_counts);

    // normalize our kmer counts
    normalize_matrix(count_matrix, 1, width);

    // multiply our kmers frequency by lambda
    for(z = 0; z < width; z++) 
      count_matrix[z] = count_matrix[z] * lambda;

    double *sensing_matrix_copy = malloc(sizeof(double) * sensing_matrix->sequences * width);
    if(sensing_matrix_copy == NULL) {
      fprintf(stderr, "Could not allocate enough memory\n");
      exit(EXIT_FAILURE);
    }

    memcpy(sensing_matrix_copy, sensing_matrix->matrix, sensing_matrix->sequences * width * sizeof(double));


    // run nnls
    double *solution = nnls(sensing_matrix_copy, count_matrix, sensing_matrix->sequences, width);

    // normalize our solution
    normalize_matrix(solution, 1, sensing_matrix->sequences);

    // add the current solution to the solutions array
    for(z = 0; z < sensing_matrix->sequences; z++ )  {
      solutions[i*sensing_matrix->sequences + z] = solution[z]; 
    }

    done++;
    printf("%ld/%ld samples processed\n", done, dir_count); 
    free(solution);
    free(count_matrix);
    free(filename);
    free(sensing_matrix_copy);
  }

  // output our matrix
  FILE *output_fh = fopen(output_filename, "w");
  if(output_fh == NULL) { 
    fprintf(stderr, "Could not open %s for writing\n", output_filename);
    exit(EXIT_FAILURE);
  }

  fprintf(output_fh, "# QIIME vQuikr OTU table\n");
  fprintf(output_fh, "#OTU_ID\t");

  // print our filename headers
  for(i = 0; i < dir_count - 1; i++) {
    fprintf(output_fh, "%s\t", filenames[i]);
  }
  fprintf(output_fh, "%s\n", filenames[dir_count - 1]);

  // get our actual values 
  for(y = 0; y < sensing_matrix->sequences; y++) {
    for(i = 0; i < dir_count; i++) {
      solutions[sensing_matrix->sequences*i + y] = round(solutions[sensing_matrix->sequences*i + y] * file_sequence_count[i]);
    }
  }

  for(y = 0; y < sensing_matrix->sequences; y++) {

    double column_sum = 0.;
    for(i = 0; i < dir_count; i++) {
      column_sum += solutions[sensing_matrix->sequences*i + y];
    }

    // if our column is zero, don't bother printing the row
    if(column_sum != 0) {
      fprintf(output_fh, "%s\t", sensing_matrix->headers[y]);

      for(i = 0; i < dir_count - 1; i++) {
				fprintf(output_fh, "%d\t", (int)solutions[sensing_matrix->sequences*i + y]);
      }
      fprintf(output_fh, "%d\n", (int)solutions[sensing_matrix->sequences*(dir_count - 1) + y]);
    }
  }
  fclose(output_fh);

  return EXIT_SUCCESS;
}
