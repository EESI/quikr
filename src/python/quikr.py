#!/usr/bin/python
import os
import sys
from StringIO import StringIO
import scipy.optimize.nnls
import scipy.sparse
import numpy as np
from subprocess import *
import gzip
import itertools

def generate_kmers(kmer):
  """ generate all possible kmers permutations seperated by newlines 

 >>> kmers =  generate_kmers(1)
 >>> generate_kmers(2)

 param kmer: the desired Mer size
 type  kmer: int
 return: Returns a string of kmers seperated by newlines
 rtype: string
 """

  return '\n'.join(''.join(x) for x in itertools.product('acgt', repeat=kmer))

def isCompressed(filename):
  """ This function checks to see if the file is gzipped
  
  >>> boolean_value = isCompressed("/path/to/compressed/gzip/file")
  >>> print boolean_value
  True

  param filename: the filename to check
  type  filename: string
  return: Returns whether the file is gzipped
  rtype: boolean

  """
  try:
    f = open(filename, "rb")
  except IOError:
    print "Warning: isCompressed could not find " + filename
    return False

  # The first two bytes of a gzipped file are always '1f 8b'
  if f.read(2) == '\x1f\x8b':
    f.close()
    return True
  else:
    f.close()
    return False

  
def train_matrix(input_file_location, kmer):
  """
  Takes a input fasta file, and kmer, returns a custom trained matrix
  """

  kmer_file_name = str(kmer) + "mers.txt"

  if not os.path.isfile(kmer_file_name):
    print "could not find kmer file"
    exit()

  input_file = Popen(["bash", "-c", "probabilities-by-read " + str(kmer) + " " + input_file_location + " <(generate_kmers 6)"], stdout=PIPE) 

  # load and  normalize the matrix by dividing each element by the sum of it's column.
  # also do some fancy rotations so that it works properly with quikr
  matrix  = np.loadtxt(input_file.stdout)
  
  matrix = np.rot90(matrix)
  matrix = matrix / matrix.sum(0)
  matrix = np.flipud(matrix);
  return matrix


def load_trained_matrix_from_file(trained_matrix_location):
  """ This is a helper function to load our trained matrix and run quikr """
  
  if isCompressed(trained_matrix_location):
    trained_matrix_file = gzip.open(trained_matrix_location, "rb")
  else:
    trained_matrix_file = open(trained_matrix_location, "rb")
  
  trained_matrix = np.load(trained_matrix_file)

  return trained_matrix


def calculate_estimated_frequencies(input_fasta_location, trained_matrix, kmer, default_lambda):
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
  # We use the count program to count 
  count_input = Popen(["count-kmers", "-r", str(kmer), "-1", "-u", input_fasta_location], stdout=PIPE) 

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
