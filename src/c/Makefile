VERSION=\"v1.0.4\"
UNAME := $(shell uname)
PWD = $(shell pwd)
CC = gcc
MULTIFASTA_CFLAGS = -pthread -L../ -I../ -std=gnu99 -fopenmp -DOMP=1
CFLAGS = -Wall -Wextra -lm -lz  -D$(UNAME) -DVERSION=$(VERSION) 


ifndef DEBUG
CFLAGS += -O3 -s -mtune=native 
else 
CFLAGS += -ggdb3 -O0 
endif

all: nnls.o kmer_utils.o quikr_functions.o quikr_train quikr multifasta_to_otu test

nnls.o: nnls.c
	$(CC) -c nnls.c -o nnls.o  $(CFLAGS)
kmer_utils.o: kmer_utils.c  quikr_functions.o
	$(CC) -c kmer_utils.c  quikr_functions.o -o kmer_utils.o  $(CFLAGS)
quikr_functions.o: quikr_functions.c 
	$(CC) -c quikr_functions.c -o quikr_functions.o  $(CFLAGS)
multifasta_to_otu: kmer_utils.o nnls.o quikr_functions.o multifasta_to_otu.c
	$(CC) multifasta_to_otu.c quikr_functions.o nnls.o kmer_utils.o -o multifasta_to_otu $(CFLAGS) $(MULTIFASTA_CFLAGS)
quikr_train: kmer_utils.o quikr_functions.o quikr_train.c
	$(CC) quikr_train.c quikr_functions.o kmer_utils.o -o quikr_train $(CFLAGS) $(QUIKR_TRAIN_CFLAGS)
quikr: kmer_utils.o nnls.o quikr_functions.o quikr.c
	$(CC) quikr.c quikr_functions.o nnls.o kmer_utils.o -o quikr $(CFLAGS) $(QUIKR_CFLAGS)
clean:
	rm -v quikr_train quikr multifasta_to_otu *.o
test: kmer_utils.o nnls.o quikr_functions.o test.c
	$(CC) test.c quikr_functions.o nnls.o kmer_utils.o -o test $(CFLAGS) -I$(PWD)
