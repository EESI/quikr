import os
import sys
import scipy.optimize.nnls
import scipy.sparse
import numpy as np
import scipy as sp
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from subprocess import *
import argparse
import platform

def main(args):

    parser = argparse.ArgumentParser(description="An argparse example")
    parser.add_argument("--trained", help="The trained matrix")
    parser.add_argument("--fasta", help="input fasta file")

    args = parser.parse_args()

    if args.trained== "install":
        print "You asked for installation" 
    else:
        print "You asked for something other than installation"

  
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
    count_input = Popen(["count-linux", "-r " + kmer, "-1", "-u", input_fasta], stdout=PIPE) 
  elif uname == "Darwin" and os.path.isfile("./count-osx"):
    print "Detected Mac OS X" 
    count_input = Popen(["count-osx", "-r 6", "-1", "-u", input_fasta], stdout=PIPE) 

  
  counts = np.loadtxt(count_input.stdout) # create a ndarray
  counts = counts / np.sum(counts) # form a probablility vector from our counts
 
  counts = default_lambda * counts

  trained_matrix  = np.loadtxt(trained_matrix_location)
   
  xstar = nnls(trained_matrix, counts)
  return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
