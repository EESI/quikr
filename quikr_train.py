#!/usr/bin/python
import numpy as np
import os
import sys
from subprocess import *
import platform
import argparse

def main():
  """
  You can call this script independently, and will save the 
  trained matrix as a numpy file.
  example: python quikr-train.py input.fasta 6 trained_matrix.npy
  """
  parser = argparse.ArgumentParser(description=
  " quikr_train returns a custom trained matrix that can be used with \
    the quikr function. \n You must supply a kmer. \n ")

  parser.add_argument("-i", "--input", help="path to the database", required=True)
  parser.add_argument("-o", "--output", help="path to output", required=True)
  parser.add_argument("-k", "--kmer", type=int, help="specifies which kmer to use", required=True)

  args = parser.parse_args()

  if not os.path.isfile(args.input):
    parser.error( "Input database not found")

  # call the quikr train function, save the output with np.save
  matrix = quikr_train(args.input, args.kmer)
  matrix = np.rot90(matrix)

  np.save(args.output, matrix)
  np.savetxt("python.csv", matrix, delimiter=",")

  return 0

def quikr_train(input_file_location, kmer):
  """
  Takes a input fasta file, and kmer, returns a custom trained matrix
  """

  kmer_file_name = str(kmer) + "mers.txt"

  if not os.path.isfile(kmer_file_name):
    print "could not find kmer file" 
    exit()

  
  uname = platform.uname()[0]

  if uname == "Linux": 
    print "Detected Linux"
    input_file = Popen(["./probabilities-by-read-linux", str(kmer), input_file_location, kmer_file_name], stdout=PIPE) 
  elif uname == "Darwin":
    print "Detected Mac OS X" 
    input_file = Popen(["./probabilities-by-read-osx", str(kmer), input_file_location, kmer_file_name]) 

  # load and  normalize the matrix by dividing each element by the sum of it's column.
  matrix  = np.loadtxt(input_file.stdout)
  
  matrix = np.rot90(matrix)
  normalized = matrix / matrix.sum(0)
  return normalized
  
if __name__ == "__main__":
    sys.exit(main())
