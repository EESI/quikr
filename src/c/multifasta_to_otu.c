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
#define USAGE "Usage:\n\tmultifasta_to_otu [OPTION...] - create a QIIME OTU table based on Quikr results. \n\nOptions:\n\n-i, --input-directory\n\tthe directory containing the samples' fasta files of reads (note each file should correspond to a separate sample)\n\n-s, --sensing-matrix\n\t location of the sensing matrix. (sensing from quikr_train)\n\n-k, --kmer\n\tspecify what size of kmer to use. (default value is 6)\n\n-l, --lambda\n\tlambda value to use. (default value is 10000)\n\n-j, --jobs\n\t specifies how many jobs to run at once. (default value is the number of CPUs)\n\n-o, --output\n\tthe OTU table, with NUM_READS_PRESENT for each sample which is compatible with QIIME's convert_biom.py (or a sequence table if not OTU's)\n\n-v, --verbose\n\tverbose mode.\n\n-V, --version\n\tprint version."

static int cmp (const void * a, const void * b) {
	return ( *(double*)a - *(double*)b );
}

char **get_fasta_files_from_file(char *fn) {
	char **files;
	int files_count = 0;
	
	// getline stuff
	ssize_t read;
	size_t len = 0;
	char *line = NULL;

  FILE *fh = fopen(fn, "r");
  if(fh == NULL) {
    fprintf(stderr, "Error opening %s - %s\n", fn, strerror(errno));
    exit(EXIT_FAILURE);
  }

	files = malloc(sizeof(char **));

	while ((read = getline(&line, &len, fh)) != -1) {
		char *file = malloc(sizeof(char) * (strlen(line)));	
		if(file == NULL) {
			exit(EXIT_FAILURE);
		}
		strncpy(file, line, strlen(line) + 1 );
		printf("%zu\n", strlen(line));
		printf("%zu\n", strlen(file));
		file[strlen(file)- 1] = '\0';
		printf(">%s< %zu >%s< %zu\n", line, strlen(line), file, strlen(file));
	
		if(access(file, F_OK) == 0) {

			files[files_count] = file;
			files_count++;

			files = realloc(files, sizeof(char *) * (files_count + 1));
			if(files == NULL) {
				fprintf(stderr, "could not realloc keys\n");
				exit(EXIT_FAILURE);
			}
		} 
		else {
			fprintf(stderr, "Warning: ignoring %s (%s)\n", file, strerror(errno)); 
			errno = 0;
			free(file);
		}
	}

	files_count++;
	files = realloc(files, sizeof(char *) * files_count);
	if(files == NULL) {
		fprintf(stderr, "could not realloc keys\n");
		exit(EXIT_FAILURE);
	}
	files[files_count] = NULL;
	
	fclose(fh);
	return files;
}

char **get_fasta_files_from_directory(char *directory) {

	DIR *dh;
	struct dirent *e;
	char **headers;
	long long count = -2; // -2 for ../ and ./
	long long i = 0;


	// open our directory
	dh = opendir(directory);
	if(dh == NULL) {
		fprintf(stderr, "could not open %s\n", directory);
		exit(EXIT_FAILURE);
	}

	while((e = readdir(dh))) 
		count++;

	e = NULL;
	rewinddir(dh);

	if(count == 0) {
		fprintf(stderr, "%s is empty\n", directory);
		exit(EXIT_FAILURE);
	}

	headers = malloc(count * sizeof(char *));
	check_malloc(headers, NULL);


	int array_pos = 0;
	for(i = 0; i < count; i++) {
		char *ext = NULL;
		e = readdir(dh);

		if(strcmp(e->d_name, "..") == 0 || strcmp(e->d_name, ".") == 0) {
			i--;
			continue;
		}

		ext = strrchr(e->d_name, '.');

		if(str_eq(ext, ".fasta") || 
				str_eq(ext, ".fa") ||
				str_eq(ext, ".fna"))
		{

			char *header = malloc(strlen(directory) + strlen(e->d_name) + 1);
			check_malloc(header, NULL);
			sprintf(header, "%s/%s", directory, e->d_name); 
			headers[array_pos] = header;
		}
		else {
			continue;
		}

		array_pos++;
	}

	headers[array_pos] = NULL;

	closedir(dh);
	return headers;
}


int main(int argc, char **argv) {

	int c;

	char *input_fasta_directory = NULL;
	char *input_fasta_filelist = NULL;
	char *sensing_matrix_filename = NULL;
	char *output_filename = NULL;

	unsigned long long i = 0;
	unsigned long long j = 0;

	unsigned long long width = 0;

	unsigned int kmer = 6;
	unsigned long long lambda = 10000;

	double rare_percent = 1.0;

	unsigned int jobs = 1;
	long done = 0;

	unsigned long long dir_count = 0;

	#ifdef Linux
		jobs = get_nprocs();
	#endif

	#ifdef Darwin
		jobs = sysconf (_SC_NPROCESSORS_ONLN);
	#endif

	int verbose = 0;

	static struct option long_options[] = {
		{"input-directory", required_argument, 0, 'i'},
		{"input-filelist", required_argument, 0, 'f'},
		{"kmer",  required_argument, 0, 'k'},
		{"lambda",  required_argument, 0, 'l'},
		{"jobs",  required_argument, 0, 'j'},
		{"output", required_argument, 0, 'o'},
		{"sensing-matrix", required_argument, 0, 's'},
		{"rare-percent", required_argument, 0, 'r'},
		{"verbose", no_argument, 0, 'v'},
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'V'},
		{0, 0, 0, 0}
	};

	while (1) {
		int option_index = 0;

		c = getopt_long (argc, argv, "f:k:l:s:i:o:j:r:hvV", long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
			case 'i':
				input_fasta_directory = optarg;
				break;
			case 'f':
				input_fasta_filelist = optarg;
				break;
			case 'j':
				jobs = atoi(optarg);
				break;
			case 'k':
				kmer = atoi(optarg);
				break;
			case 'l':
				lambda = atoi(optarg);
				break;
			case 'o':
				output_filename = optarg;
				break;
			case 'r':
				rare_percent = atof(optarg);
				break;
			case 's':
				sensing_matrix_filename = optarg;
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

 	// input fasta parsing
	if(input_fasta_directory == NULL && input_fasta_filelist == NULL) {
		fprintf(stderr, "Error: input fasta directory (-i) or input fasta filelist (-f) must be specified\n\n");
		fprintf(stderr, "%s\n", USAGE);
		exit(EXIT_FAILURE);
	}

	if(input_fasta_directory != NULL && input_fasta_filelist != NULL) {
		fprintf(stderr, "Error: input fasta directory (-i) and input fasta filelist (-f) cannot be used concurrently\n\n");
		fprintf(stderr, "%s\n", USAGE);
		exit(EXIT_FAILURE);
	}

	if(rare_percent <= 0 || rare_percent > 1.0) {
		fprintf(stderr, "Error: rare percent must be between 0 and 1\n");
		exit(EXIT_FAILURE);
	}

	if(verbose) { 
		printf("kmer: %u\n", kmer);
		printf("lambda: %llu\n", lambda);
		printf("input directory: %s\n", input_fasta_directory);
		printf("input filelist: %s\n", input_fasta_filelist);
		printf("sensing database: %s\n", sensing_matrix_filename);
		printf("output: %s\n", output_filename);
		printf("number of jobs to run at once: %d\n", jobs); 
	}

	if(access (sensing_matrix_filename, F_OK) == -1) {
		fprintf(stderr, "Error: could not find %s\n", sensing_matrix_filename);
		exit(EXIT_FAILURE);
	}

	if(kmer == 0) {
		fprintf(stderr, "Error: zero is not a valid kmer\n");
		exit(EXIT_FAILURE);
	}

	// load filenames
	char **filenames = NULL;
	if(input_fasta_directory != NULL)
		filenames = get_fasta_files_from_directory(input_fasta_directory);
	else
		filenames = get_fasta_files_from_file(input_fasta_filelist);

	while(filenames[dir_count] != NULL)
		dir_count++;
	
	if(dir_count == 0) {
		fprintf(stderr, "Error: No files loaded from input\n");
		exit(EXIT_FAILURE);
	}

	// 4 "ACGT" ^ Kmer gives us the size of output rows
	width = pow(4, kmer);

	struct matrix *sensing_matrix = load_sensing_matrix(sensing_matrix_filename, kmer);
	double *sensing_matrix_ptr = sensing_matrix->matrix;
	unsigned long long sequences = sensing_matrix->sequences;

	if(verbose) {
		printf("directory count: %llu\n", dir_count);
		printf("width: %llu\n", width);
		printf("sequences: %llu\n", sequences);
	}

  unsigned long long *solutions = malloc(dir_count * sequences * sizeof(unsigned long long));
	check_malloc(solutions, NULL);

  long long *file_sequence_count = calloc(dir_count, sizeof(long long));
	check_malloc(file_sequence_count, NULL);

	#ifdef OMP
  	omp_set_num_threads(jobs);
	#endif


  printf("Beginning to process samples\n"); 

  #pragma omp parallel for shared(solutions, sensing_matrix_ptr, width, done, sequences)
  for(size_t i = 0; i < dir_count; i++ ) {

		size_t x = 0;
		size_t y = 0;
		size_t z = 0;

		unsigned long long file_sequence_count = 0;
  	unsigned long long rare_value = 0;
  	unsigned long long rare_width = 0; 

		double rare_percent = 1.0;
		
		printf("processing %s\n", filenames[i]);
		file_sequence_count = count_sequences(filenames[i]);

		// load counts matrix
		double *count_matrix = malloc(width * sizeof(double));
		check_malloc(count_matrix, NULL);
		double *sorted_count_matrix = malloc(width * sizeof(double));
		check_malloc(sorted_count_matrix, NULL);

		// convert our matrix into  doubles
		{
			unsigned long long *integer_counts = get_kmer_counts_from_file(filenames[i], kmer);

			for(x = 0; x < width; x++) {
				sorted_count_matrix[x] = count_matrix[x] = (double)integer_counts[x];
			}

			free(integer_counts);
		}

		// sort our array
		qsort(sorted_count_matrix, width, sizeof(double), cmp);

		// get our "rare" counts
		for(y = 0; y < width; y++) {
			double percentage = 0;

			rare_value = sorted_count_matrix[y];
			rare_width = 0;
			for(x = 0; x < width; x++) {
				if(count_matrix[x] <= rare_value) {
					rare_width++;
				}
			}
			percentage = (double)rare_width / (double)width;

			if(percentage >= rare_percent)
				break;
		}
		
		free(sorted_count_matrix);

		if(verbose)
			printf("there are %llu values less than %llu\n", rare_width, rare_value);

		// add a extra space for our zero's array
		rare_width++;

		// store our count matrix
		double *count_matrix_rare = calloc(rare_width, sizeof(double));
		check_malloc(count_matrix_rare, NULL);

		double *sensing_matrix_rare = malloc((rare_width) * sequences * sizeof(double));
		check_malloc(sensing_matrix_rare, NULL);

		// copy only kmers from our original counts that match our rareness percentage
		// in both our count matrix and our sensing matrix
		for(x = 0, y = 1;  x < width; x++) {
			if(count_matrix[x] < rare_value) {
				count_matrix_rare[y] = count_matrix[x];

				for(z = 0; z < sequences; z++) 
					sensing_matrix_rare[z*rare_width + y] = sensing_matrix_ptr[z*width + x];

				y++;
			}
		}

		// normalize our kmer counts and our sensing_matrix
		normalize_matrix(count_matrix_rare, 1, rare_width);
		normalize_matrix(sensing_matrix_rare, sequences, rare_width);

		// multiply our kmer counts and sensing matrix by lambda
		for(x = 1; x < rare_width; x++) 
			count_matrix_rare[x] *= lambda;

		for(x = 0; x < sequences; x++) {
			for(y = 1; y < rare_width; y++) {
				sensing_matrix_rare[rare_width*x + y] *= lambda;
			}
		}

		// count_matrix's first element should be zero
		count_matrix_rare[0] = 0;
		// stack one's on our first row of our sensing matrix
		for(x = 0; x < sequences; x++) {
			sensing_matrix_rare[x*rare_width] = 1.0;
		}

		double *solution = nnls(sensing_matrix_rare, count_matrix_rare, sequences, rare_width);

		// add the current solution to the solutions array
		for(unsigned long long z = 0; z < sequences; z++ )  {
			solutions[sensing_matrix->sequences*i + z] = (unsigned long long)round(solution[z] * file_sequence_count);
		}

		done++;
		printf("%ld/%llu samples processed\n", done, dir_count); 
		free(solution);
		free(count_matrix_rare);
		free(count_matrix);
		free(sensing_matrix_rare);
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

	for(j = 0; j < sequences; j++) {

		double column_sum = 0.;
		for(i = 0; i < dir_count; i++) {
			column_sum += solutions[sensing_matrix->sequences*i + j];
		}

		// if our column is zero, don't bother printing the row
		if(column_sum != 0) {
			fprintf(output_fh, "%s\t", sensing_matrix->headers[j]);

			for(i = 0; i < dir_count - 1; i++) {
				fprintf(output_fh, "%llu\t", solutions[sensing_matrix->sequences*i + j]);
			}
			fprintf(output_fh, "%llu\n", solutions[sensing_matrix->sequences*(dir_count - 1) + j]);
		}
	}
	fclose(output_fh);

	return EXIT_SUCCESS;
}
