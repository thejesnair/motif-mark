#!/usr/bin/env python

'''
Docstring for motif-mark.motif-mark-oop
Takes in FASTA file and generates image marking motifs along gene
Script outputs both an svg and png
'''

from bioinfo import oneline_fasta # converts fasta record to oneline
from bioinfo import validate_base_seq # determine if RNA (true) or DNA seq

import argparse
import re # regex

def get_args():
    '''Takes in command line arguments for paths to motifs files, FASTA file, and output image.'''
    parser = argparse.ArgumentParser(description=
                                     "Program to generate an image of motifs along a gene")
    #parser.add_argument("-f", "--file", help="path for fasta file", type=str, required=True)
    parser.add_argument("-m", "--motif", help="deduped SAM file, absolute path", type=str, required=True)
    #parser.add_argument("-o", "--outfile", help = "name for outputted motif image (png)", type=str, required=True)

    return parser.parse_args()

args = get_args()
#fasta = args.file
motif_file = args.motif
#image = args.outfile


# motifs set
with open(motif_file) as f:
    motifs = {line.strip() for line in f} # strips trailing space, adds motifs to SET

print(motifs) #check that properly pulls motifs
