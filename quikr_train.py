#from scipy.sparse import * 
import numpy as np
import sys
from subprocess import *
import platform

# You can call this script independently, and will save the 
# trained matrix as a numpy file.
# example: python quikr-train.py input.fasta 6 trained_matrix.npy

def main(argv):
  input_file_location = argv[1]
  kmer = argv[2]
  output_file_location = argv[3]

  # call the quikr train function, save the output with np.save
  matrix = quikr_train(argv[1], argv[2])
  np.save(output_file_location, matrix)

  return 0

def quikr_train(input_file_location, kmer):

  
  print "input fasta training file: " + input_file_location
  print "kmer: " + kmer

  kmer_file_name = kmer + "mers.txt"
  print kmer_file_name

  
  uname = platform.uname()[0]

  if uname == "Linux": 
    print "Detected Linux"
    input_file = Popen(["./probabilities-by-read-linux", kmer, input_file_location, kmer_file_name], stdout=PIPE) 
  elif uname == "Darwin":
    print "Detected Mac OS X" 
    input_file = Popen(["./probabilities-by-read-osx", kmer, input_file_location, kmer_file_name]) 

  # load and  normalize the matrix by dividing each element by the sum of it's column.
  matrix  = np.loadtxt(input_file.stdout)
  normalized = matrix / matrix.sum(0)

  return normalized
  


if __name__ == "__main__":
    sys.exit(main(sys.argv))
