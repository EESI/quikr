#!/usr/bin/python
from multiprocessing import Pool
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
trained_matrix = ""
output_directory = ""


def main():


    #do: write up the description
    parser = argparse.ArgumentParser(description="MultifastaOTU")

    parser.add_argument("-i", "--input", help="directory containing fasta files", required=True)
    parser.add_argument("-o", "--otu-table", help="otu_table", required=True)
    parser.add_argument("-t", "--trained-matrix", help="otu_table", required=True)
    parser.add_argument("-d", "--output-directory", help="quikr output directory", required=True)
    parser.add_argument("-l", "--lamb", type=int, help="the default lambda value is 10,000")
    parser.add_argument("-k", "--kmer", type=int, help="specifies which kmer to use, default=6")
    parser.add_argument("-j", "--jobs", type=int, help="specifies how many jobs to run at once, default=number of CPUs")
    args = parser.parse_args()
    
    # our defaults
    trained_matrix = args.trained_matrix

    # Make sure our input exist
    if not os.path.isdir(args.input):
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
    fasta_list = os.listdir(args.
    pool = Pool(processes=jobs)
    result = pool.map(quikr_call, fasta_list)
    return 0

def quikr_call(fasta_file):
  xstar = q.quikr(fasta_file, training_matrix, kmer, lamb)
  np.savetxt(output_directory + os.path.basename(fasta_file), xstar, delimiter=",", fmt="%f")
  return 0

 if __name__ == "__main__":
     sys.exit(main())
 
