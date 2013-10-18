VERSION=\"v1.0.4\"
UNAME := $(shell uname)
CC = gcc
QUIKR_TRAIN_CFLAGS = -D$(UNAME) -DVERSION=$(VERSION)
MULTIFASTA_CFLAGS = -pthread -L../ -I../ -std=gnu99 -fopenmp -D$(UNAME) -DVERSION=$(VERSION)
QUIKR_CFLAGS = -std=gnu99 -D$(UNAME) -DVERSION=$(VERSION)

ifndef DEBUG
CFLAGS = -O3 -s -mtune=native -Wextra -Wall -lm -lz
else 
CFLAGS = -ggdb3 -O0 -Wall -Wextra -lm -lz
endif

all: nnls.o kmer_utils.o quikr_train quikr multifasta_to_otu

nnls.o: nnls.c
	$(CC) -c nnls.c -o nnls.o  $(CFLAGS)
kmer_utils.o: kmer_utils.c 
	$(CC) -c kmer_utils.c -o kmer_utils.o  $(CFLAGS)
multifasta_to_otu: kmer_utils.o nnls.o multifasta_to_otu.c
	$(CC) multifasta_to_otu.c quikr_functions.c nnls.o kmer_utils.o -o multifasta_to_otu $(CFLAGS) $(MULTIFASTA_CFLAGS)
quikr_train: kmer_utils.o quikr_train.c
	$(CC) quikr_train.c quikr_functions.c -o quikr_train $(CFLAGS) $(QUIKR_TRAIN_CFLAGS)
quikr: kmer_utils.o nnls.o quikr.c
	$(CC) quikr.c quikr_functions.c nnls.o kmer_utils.o -o quikr $(CFLAGS) $(QUIKR_CFLAGS) -I../
clean:
	rm -vf quikr_train quikr multifasta_to_otu *.o
