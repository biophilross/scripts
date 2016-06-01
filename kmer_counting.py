#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Description: Read in fasta file line by line and record every kmer
as defined by the user
"""

# import modules
from __future__ import division
import sys
import argparse
import operator

def main():

    parser = argparse.ArgumentParser( \
                                     description = "Compute metrics of kmers within fasta file")
    parser.add_argument("-f", dest = "fasta_file", help = "Fasta file")
    parser.add_argument("-k", dest = "kmer_length", help = "Kmer length")
    if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)
    args = parser.parse_args()

    fasta_file = args.fasta_file
    kmer_len = args.kmer_length

    #---- STEP1: build a hash with all kmers from the reads ----#
    seqs = 0
    kmers = dict()

    with open(fasta_file) as fasta:
        for line in fasta:
            if(line[0] == ">"):
                seqs += 1
                pass
            else:
                seq = line.rstrip()
                for i in range(len(seq)):
                    if(int(kmer_len) + i <= len(seq)):
                        kmer = seq[i:int(kmer_len) + i]
                        if(kmer not in kmers):
                            kmers[kmer] = 1
                        elif(kmer in kmers):
                            kmers[kmer] += 1

    kmers = sorted(kmers.iteritems(), key=operator.itemgetter(1), reverse=True)

    count = 0
    running_total = 0
    top_5 = ""
    bottom_5 = ""
    for kmer in kmers:
        count += 1
        running_total += kmer[1]
        if(count <= 5):
            top_5 += "%s\t%s\n" % (kmer[0], kmer[1])
        if(count <= 10):
            at_10 = running_total
        if(count <= 100):
            at_100 = running_total
        if(count <= 1000):
            at_1000 = running_total
        if(count <= 10000):
            at_10000 = running_total
    for i in range(1, 6):
        least_common_kmer = kmers[-i]
        bottom_5 += "%s\t%s\n" % (least_common_kmer[1], least_common_kmer[0])

    #---- STEP3: report the metrics ----#
    print >> sys.stdout, ""
    print >> sys.stdout, "Used %s to find kmer size %s metrics" % ("program", int(kmer_len))
    print >> sys.stdout, "Reads file: %s\n" % ("reads")
    print >> sys.stdout, "Total number of reads: %s" % (reads)
    print >> sys.stdout, "Total number of kmers: %s\n" % (len(kmers) + 5)
    print >> sys.stdout, "Top 10 kmers combined:    %s     %.2f percent of all kmers" % (str(at_10) + "/" + str(running_total), (at_10/running_total)*100)
    print >> sys.stdout, "Top 100 kmers combined:   %s     %.2f percent of all kmers" % (str(at_100) + "/" + str(running_total), (at_100/running_total)*100)
    print >> sys.stdout, "Top 1000 kmers combined:  %s     %.2f precent of all kmers" % (str(at_1000) + "/" + str(running_total), (at_1000/running_total)*100)
    print >> sys.stdout, "Top 10000 kmers combined: %s     %.2f precent of all kmers\n" % (str(at_10000) + "/" + str(running_total), (at_10000/running_total)*100)
    print >> sys.stdout, "5 most common kmers:\n%s\n5 least common kmers:\n%s\n" % (top_5, bottom_5)

if __name__ == "__main__":
    main()
