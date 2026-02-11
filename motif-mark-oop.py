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


# motifs list
with open(motif_file) as f:
    motifs = [line.strip() for line in f] # strips trailing space, adds motifs to SET

print(motifs) #check that properly pulls motifs

# create Motif objects 
class Motif:
    ''' Holds one motif definition, regex pattern, ambiguity rules 
        What is the motif I'm looking for? '''

    # degenerate bases (ambiguity)
    ''' Y: pyrimidine
        R: purine
        N: any
        V: not T (or U)
        H: not G
        D: not C
        B: not A
        Note: for RNA, replace T with U'''

    DNA_degenerate = {'Y': "CT", 'R': "AG", 'N': "ACGT",
                    'V': "ACG", 'H': 'ACT', 'D': "AGT",
                    'B': "CGT"}
    
    RNA_degenerate = {'Y': "CU", 'R': "AG", 'N': "ACGU",
                    'V': "ACG", 'H': 'ACU', 'D': "AGU",
                    'B': "CGU"}
    
    def __init__(self, motif, mode):

        ## DATA
        self.motif = motif
        self.mode = mode

        # normalize
        self.pattern = motif.upper() # normalize for checking
        if mode == "DNA":
            self.degenerate_map = self.DNA_degenerate # have to reference class variables
            if "U" in self.pattern:
                raise ValueError("Motif contains U but sequence is DNA")
        elif mode == "RNA":
            self.degenerate_map = self.RNA_degenerate
            if "T" in self.pattern:
                raise ValueError("Motif contains T but sequence is RNA")
        else:
            raise ValueError("Input must be RNA or DNA")

        self.length = len(self.pattern)

        # build regex strings
        self.regex_body = None # WILL BUILD LATER
        self.lookahead_overlap = None # WILL BUILD LATER

    ## METHODS
    def build_regex_body(self):
        pass

    def build_lookahead_overlap_regex(self):
        pass


class SplicingRegion:
    pass

class MotifLocation:
    pass

class MotifScanner:
    pass

class MotifMarkRenderer:
    pass

