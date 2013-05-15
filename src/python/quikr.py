'''
  This is a python implementation of the Quikr program.
'''
__author__ = "Calvin Morrison"
__credits__ = ["Calvin Morrison", "Gail Rosen", "Simon Foucart", "Jean-Luc Bouchot", "David Koslicki"]
__license__ = "GPL"
__maintainer__ = "Gail Rosen"
__email__ = "gailro@gmail.com"
__status__ = "Release"

import os
import sys
from StringIO import StringIO
import scipy.optimize.nnls
import scipy.sparse
import numpy as np
from subprocess import *
import gzip
import itertools


def is_compressed(filename):
  """ This is a helper function to see if a file is gzipped
  
  >>> boolean_value = is_compressed("/path/to/compressed/gzip/file")
  >>> print boolean_value
  True

  input: filename
  output: boolean
  """
  try:
    f = open(filename, "rb")
  except IOError:
    print "Warning: is_compressed could not find " + filename
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
  Takes a input fasta file, and kmer, returns a custom sensing matrix
  
  returns an ndarray
  """

  input_file = Popen(["bash", "-c", "probabilities-by-read " + str(kmer) + " " + input_file_location + " <(generate_kmers 6)"], stdout=PIPE) 

  # load and  normalize the matrix by dividing each element by the sum of it's column.
  # also do some fancy rotations so that it works properly with quikr
  matrix  = np.loadtxt(input_file.stdout)
  
  matrix = np.rot90(matrix)
  matrix = matrix / matrix.sum(0)
  matrix = np.flipud(matrix);
  return matrix


def load_sensing_matrix_from_file(sensing_matrix_location):
  """ This is a helper function to load our sensing matrix and run quikr """
  
  if is_compressed(sensing_matrix_location):
    sensing_matrix_file = gzip.open(sensing_matrix_location, "rb")
  else:
    sensing_matrix_file = open(sensing_matrix_location, "rb")
  
  sensing_matrix = np.load(sensing_matrix_file)

  return sensing_matrix


def calculate_estimated_frequencies(input_fasta_location, sensing_matrix, kmer, default_lambda):
  """
  input_fasta is the input fasta file to find the estimated frequencies of
  sensing_matrix is the sensing matrix we are using to estimate the species
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
  sensing_matrix = sensing_matrix * default_lambda;
  sensing_matrix = np.vstack((np.ones(sensing_matrix.shape[1]), sensing_matrix))

  xstar, rnorm = scipy.optimize.nnls(sensing_matrix, counts) 
  xstar = xstar / xstar.sum(0) 

  return xstar
