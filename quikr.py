#!/usr/bin/python
import os
import sys
import scipy.optimize.nnls
import scipy.sparse
import numpy as np
from subprocess import *
import argparse
import platform

def main():

    parser = argparse.ArgumentParser(description=
    "Quikr returns the estimated frequencies of batcteria present when given a \
    input FASTA file. \n \
    A default trained matrix will be used if none is supplied \n \
    You must supply a kmer and default lambda if using a custom trained \
    matrix.")

    parser.add_argument("-f", "--fasta", help="path to a fasta file", required=True)
    parser.add_argument("-t", "--trained-matrix", help="path to a custom trained matrix")
    parser.add_argument("-l", "--lamb", type=int, help="the default lambda value is 10,000")
    parser.add_argument("-k", "--kmer", type=int, 
        help="specifies which kmer to use, must be used with a custom trained database")


    args = parser.parse_args()

    # Do some basic sanity checks

    if not os.path.isfile(args.fasta):
        parser.error( "Input fasta file not found")
    
    # If we are using a custom trained matrix, we need to do some basic checks
    if args.trained_matrix is not None:  
         
        if not os.path.isfile(args.trained_matrix):
            parser.error("custom trained matrix not be found")

        if args.kmer is None:
            parser.error("A kmer is required when using a custom matrix")
        else:
          kmer = args.kmer

        if args.lamb is None:
            # use 10,000 as default Lambda
            input_lambda = 10000
    # If we aren't using a custom trained matrix, load in the defaults
    else:
        trained_matrix_location = "output.npy"
        input_lambda = 10000
        kmer = 6
        xstar = quikr(args.fasta, trained_matrix_location, kmer, input_lambda)
        
    return 0

def quikr(input_fasta_location, trained_matrix_location, kmer, default_lambda):
  """
  input_fasta is the input fasta file to find the estimated frequencies of
  trained_matrix is the trained matrix we are using to estimate the species
  kmer is the desired k-mer to use
  default_lambda is inp 
  
  returns the estimated requencies of bacteria present when given an input
  FASTA file of amplicon (454) reads. A k-mer based, L1 regularized, sparsity
  promoting algorthim is utilized. 

  In practice reconstruction is accurate only down to the genus level (not 
  species or strain).
  """

  uname = platform.uname()[0]

  # We use the count program to count ____
  if uname == "Linux" and os.path.isfile("./count-linux"):
    print "Detected Linux"
    count_input = Popen(["./count-linux", "-r", str(kmer), "-1", "-u", input_fasta_location], stdout=PIPE) 
  elif uname == "Darwin" and os.path.isfile("./count-osx"):
    print "Detected Mac OS X" 
    count_input = Popen(["count-osx", "-r", str(kmer), "-1", "-u", input_fasta_location], stdout=PIPE) 

  
  # load the output of our count program and form a probability vector from the counts  
  counts = np.loadtxt(count_input.stdout) 
  counts = counts / np.sum(counts) 
 
  counts = default_lambda * counts

  trained_matrix  = np.load(trained_matrix_location)
  
  # perform the non-negative least squares
  # import pdb; pdb.set_trace()
  counts = np.rot90(counts)
  xstar = scipy.optimize.nnls(trained_matrix, counts) 
 
  xstar = xstar / sum(xstar) 
  return xstar


if __name__ == "__main__":
    sys.exit(main())
