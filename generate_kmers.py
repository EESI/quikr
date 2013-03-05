#!/usr/bin/python
import itertools
import sys

kmer_array ='\n'.join(''.join(x) for x in itertools.product('acgt', repeat=int(sys.argv[1])))

print kmer_array
