#!/usr/bin/python
import itertools
import sys

print '\n'.join(''.join(x) for x in itertools.product('acgt', repeat=int(sys.argv[1])))
