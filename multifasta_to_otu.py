#!/usr/bin/python
from multiprocessing import Pool
from Bio import SeqIO
import multiprocessing
import os
import quikr_train as qt
import quikr as q
import sys
import numpy as np
import argparse
import platform

# our defaults
kmer = 6
lamb = 10000
output_directory = ""
input_directory = ""


def main():

  global kmer
  global input_directory 
  global output_directory 
  global lamb
  global trained_matrix 
  #do: write up the description
  parser = argparse.ArgumentParser(description="MultifastaOTU")

  parser.add_argument("-i", "--input-directory", help="directory containing fasta files", required=True)
  parser.add_argument("-o", "--otu-table", help="otu_table", required=True)
  parser.add_argument("-t", "--trained-matrix", help="your trained matrix ", required=True)
  parser.add_argument("-f", "--trained-fasta", help="the fasta file used to train your matrix", required=True)
  parser.add_argument("-d", "--output-directory", help="quikr output directory", required=True)
  parser.add_argument("-l", "--lamb", type=int, help="the default lambda value is 10,000")
  parser.add_argument("-k", "--kmer", type=int, help="specifies which kmer to use, default=6")
  parser.add_argument("-j", "--jobs", type=int, help="specifies how many jobs to run at once, default=number of CPUs")
  args = parser.parse_args()
    
  # our defaults
  jobs=multiprocessing.cpu_count()
  trained_matrix = args.trained_matrix
  input_directory = args.input_directory
  output_directory = args.output_directory

  # Make sure our input exist
  if not os.path.isdir(args.input_directory):
    parser.error( "Input directory not found")

  if not os.path.isdir(args.output_directory):
    parser.error( "Input directory not found")

  if not os.path.isdir(args.output_directory):
    os.path.mkdir(args,output_directory)

  if not os.path.isfile(args.trained_matrix):
    parser.error("custom trained matrix not found")
    
  # use alternative lambda
  if args.lamb is not None:
    lamb = args.lamb
    
  if args.jobs is not None:
    jobs = args.jobs

  if args.kmer is not None:
    kmer = args.kmer

  # Load trained matrix
  trained_matrix = np.load(args.trained_matrix);

  # Return a list of the input directory
  fasta_list = os.listdir(args.input_directory)

  # Queue up and run our quikr functions.
#  pool = Pool(processes=jobs)
#  results = pool.map(quikr_call, fasta_list)

  # Create an array of headers
  records = []
  trained_matrix_headers = open(args.trained_fasta, "rU")
  for record in SeqIO.parse(trained_matrix_headers, "fasta"):
    records.append(record.id)
  trained_matrix_headers.close()

  final_output = np.zeros((len(records), len(fasta_list)))
  print len(fasta_list)

  # load the keys with values from each fasta result
  for fasta, fasta_it in map(None, fasta_list, range(len(fasta_list))):
    fasta_file = open(input_directory + fasta, "rU")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    number_of_sequences = len(sequences)
    fasta_file.close()

    print number_of_sequences
    
    proportions = np.loadtxt(output_directory + fasta);
    for proportion, proportion_it in map(None, proportions, range(len(proportions))):
     if(round(proportion * number_of_sequences) is not 0):
        print str(fasta_it) + " " + str(proportion_it)
        final_output[fasta_it, proportion_it] = proportion * number_of_sequences
  
  np.savetxt(args.otu_table, final_output, delimiter=",", fmt="%d")
      



  # Write the otu table
  return 0

def quikr_call(fasta_file):
  inp = input_directory + fasta_file
  output = output_directory + os.path.basename(fasta_file)

  xstar = q.quikr(inp, trained_matrix, kmer, lamb)
  np.savetxt(output, xstar, delimiter=",", fmt="%f")
  return xstar

if __name__ == "__main__":
  sys.exit(main())
 
