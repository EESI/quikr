#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "quikr.h"
#include "quikr_functions.h"


static int failed = 0;

void status(int test, int pass, char *testname) {
	fprintf(stdout, "test %d in %s ", test, testname);

	if(pass)
		fprintf(stdout, "passed.\n");
	else {
		fprintf(stdout, "failed. \n");
		failed = 1;
	}
}

#define test_eq(data, value) status(test_number, (data == value), test_name); test_number++; 
#define test_n_eq(data, value) status(test_number, (data != value), test_name); test_number++; 
#define test_gt(data, value) status(test_number, (data < value), test_name); test_number++; 
#define test_lt(data, value) status(test_number, (data > value), test_name); test_number++; 

#define header(test_name) printf("testing %s:\n", test_name);
#define footer() printf("\n\n");;

void test_count_sequences() {

	int test_number = 1;
	char *test_name = "test_count_sequences";

	// test 1
	// Make sure we can count this properly
	int t_count_sequences_res = count_sequences("tests/data/test.fa");
	test_eq(t_count_sequences_res, 17414);

	// test 2
	// Should be zero, since the file is empty
	t_count_sequences_res = count_sequences("tests/data/empty.fa");
	test_eq(t_count_sequences_res, 0);

	// test 3
	// Doesn't exist, count should be zero
	t_count_sequences_res = count_sequences("non-existant-file");
	test_eq(t_count_sequences_res, 0);

}


void test_normalize_matrix() {

	int test_number = 1;
	int fail_flag = 0;
	char *test_name = "test_normailze_matrix";

	double sum = 0;
	int i = 0;
	int j = 0;

	double matrix_1[4] = {1,1,1,1};
	double matrix_2[8] = {2,2,4,4,8,8,16,16};


	// test 1
	// our matrix to sum to 1 properly even with only one row 
	normalize_matrix(matrix_1, 1, 4);
	for(i = 0; i < 4; i++) 
		sum += matrix_1[i];

	test_eq(sum, 1);

	// test 2
	// each element in the array should be .25
	for(i = 0; i < 4; i++)
		if(matrix_1[i] != 0.25)
			fail_flag = 1;

	test_eq(fail_flag, 0);

	fail_flag = 0;
	sum = 0;
	// test 3
	//  make sure a multidimensional array sums each row properly
	normalize_matrix(matrix_2, 2, 4);
	for(j = 0; j < 2; j++) {
		double row_sum = 0;

		for(i = 0; i < 4; i++) 
			row_sum += matrix_2[j*2 + i];	
		
		// floating point hack
		if(fabsl(row_sum - 1.0) > .0001)
			fail_flag = 1;
	}

	test_eq(fail_flag, (double)0);
		
}

int main() {

	header("count_sequences");
	test_count_sequences();
	footer();

	header("normalize_matrix");
	test_normalize_matrix();
	footer();

	if(failed)
		return EXIT_FAILURE;
	else
		return EXIT_SUCCESS;
};
