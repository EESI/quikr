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

    parser.add_argument("-f", "--fasta", help="fasta file", required=True)
    parser.add_argument("-o", "--output", help="output path (csv output)", required=True)
    parser.add_argument("-t", "--trained-matrix", help="trained matrix", required=True)
    parser.add_argument("-l", "--lamb", type=int, help="the default lambda value is 10,000")
    parser.add_argument("-k", "--kmer", type=int, required=True,
        help="specifies which kmer to use, must be used with a custom trained database")


    args = parser.parse_args()
    
    # our default lambda is 10,000
    lamb = 10000


    # Make sure our input exist
    if not os.path.isfile(args.fasta):
        parser.error( "Input fasta file not found")

    if not os.path.isfile(args.trained_matrix):
        parser.error("custom trained matrix not found")
    
    # use alternative lambda
    if args.lamb is not None:
        lamb = args.lamb
    
    xstar = quikr_load_trained_matrix_from_file(args.fasta, args.trained_matrix, args.kmer, lamb)
    np.savetxt(args.output, xstar, delimiter=",", fmt="%f")
    return 0

def quikr_load_trained_matrix_from_file(input_fasta_location, trained_matrix_location, kmer, default_lambda):
  
  trained_matrix = np.load(trained_matrix_location)
  xstar = quikr(input_fasta_location, trained_matrix, kmer, default_lambda)
  return xstar
  
def quikr(input_fasta_location, trained_matrix, kmer, default_lambda):
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
  counts = counts / counts.sum(0) 
  counts = default_lambda * counts
  counts = np.concatenate([np.zeros(1), counts])

  #form the k-mer sensing matrix
  trained_matrix = trained_matrix * default_lambda;
  trained_matrix = np.vstack((np.ones(trained_matrix.shape[1]), trained_matrix))


  xstar, rnorm = scipy.optimize.nnls(trained_matrix, counts) 

  xstar = xstar / xstar.sum(0) 

  return xstar


if __name__ == "__main__":
    sys.exit(main())
