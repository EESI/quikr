VERSION=\"v1.0.3\"
UNAME := $(shell uname)
CC = gcc
QUIKR_TRAIN_CFLAGS = -O3 -s -mtune=native -Wall -lm -lz -D$(UNAME) -DVERSION=$(VERSION)
MULTIFASTA_CFLAGS = -O3 -s -mtune=native -Wextra -Wall -lm -lz -pthread -L../ -I../ -std=gnu99 -fopenmp -D$(UNAME) -DVERSION=$(VERSION)
QUIKR_CFLAGS= -O3 -s -mtune=native -Wextra -Wall -lm -lz -std=gnu99 -D$(UNAME) -DVERSION=$(VERSION)

all: quikr_train_ quikr_ multifasta_to_otu_

multifasta_to_otu_:
	$(CC) multifasta_to_otu.c quikr_functions.c nnls.c -o multifasta_to_otu $(MULTIFASTA_CFLAGS)
quikr_train_:
	$(CC) quikr_train.c quikr_functions.c -o quikr_train $(QUIKR_TRAIN_CFLAGS)
quikr_:
	$(CC) quikr.c quikr_functions.c nnls.c -o quikr $(QUIKR_CFLAGS) -I../
clean:
	rm -vf quikr_train quikr multifasta_to_otu
