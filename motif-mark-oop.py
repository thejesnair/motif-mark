#!/usr/bin/env python

'''
Docstring for motif-mark.motif-mark-oop
Takes in FASTA file and generates image marking motifs along gene
Script outputs both an svg and png
'''

import argparse
import re # regex

def get_args():
    '''Takes in command line arguments for paths to motifs files, FASTA file, and output image.'''
    parser = argparse.ArgumentParser(description=
                                     "Program to generate an image of motifs along a gene")
    parser.add_argument("-f", "--file", help="path for fasta file", type=str, required=True)
    parser.add_argument("-m", "--motif", help="deduped SAM file, absolute path", type=str, required=True)
    #parser.add_argument("-o", "--outfile", help = "name for outputted motif image (png)", type=str, required=True)

    return parser.parse_args()

args = get_args()
input_fasta = args.file
motif_file = args.motif
#image = args.outfile

## HELPER FXN ##

# read in fasta
def read_fasta(input_fa:str) -> list[tuple[str,str]]:
    ''' Parses a FASTA file and returns a list of (header, sequence) tuples
        Handles multiline sequences and stores seq as oneline '''
    
    current_header = None
    current_seq = []
    records = []
    
    with open(input_fa, 'r') as in_fa:
        for line in in_fa:
            line = line.strip()
            if line.startswith('>'):    # new header, new record starting
                if current_header is not None:      # if there is already have a record in store, save
                    records.append((current_header, ''.join(current_seq)))

                # start tracking new record
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)        # reading seq lines for the current record
        
        # Saving last record
        if current_header is not None:
            records.append((current_header, ''.join(current_seq)))
    return records

records = read_fasta(input_fasta)

# motifs list
with open(motif_file) as f:
    motifs = [line.strip() for line in f] # strips trailing space, adds motifs to list

print(motifs) #check that properly pulls motifs

## CLASSES ##
 
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
        self.regex_body = self.build_regex_body()       # call method
        self.lookahead_overlap = self.build_lookahead_overlap_regex()

    ## METHODS
    def build_regex_body(self) -> str:
        ''' Build regex expression for motif '''
        regex_exp = []
        for char in self.pattern.upper():       # upper bc of dict
            if char in self.degenerate_map:
                regex_exp.append("[" + self.degenerate_map[char] + "]")
            else:
                regex_exp.append(char)
        return ''.join(regex_exp)

    def build_lookahead_overlap_regex(self) -> str:
        ''' Build regex expression for lookahead/overlap '''
        return "(?=(" + self.regex_body + "))"      # (?=...) is a lookahead assertion!

class SplicingRegion:
    ''' Represents FASTA record
        Where I'm looking for motifs '''
    
    def __init__(self, header: str, sequence: str):
        self.header = header
        self.sequence = sequence       # original sequence
        self.sequence_upper = sequence.upper()    # to identify motifs

        self.introns = []
        self.exons = []

        self.mode = self.detect_mode()
        self.find_intron_exon_regions()

    ## METHODS
    def detect_mode(self) -> str:
        ''' Determine if DNA or RNA '''
        if "U" in self.sequence_upper:
            return "RNA"
        elif "T" in self.sequence_upper:
            return "DNA"
        else:
            raise ValueError("DNA or RNA FASTA required")
        

    def find_intron_exon_regions(self):
        ''' Determines coordinates of introns and exons
            Tuple: Intron/exon number, start pos, end pos '''

        intron_num = 0
        exon_num = 0
        seq = self.sequence
        n = len(seq)
        i = 0       # start coordinate
        j = 0       # 'scanning' end coordinate

        while i < n:
            if seq[i].islower():
                j = i       # start j at current i location
                while j < n and seq[j].islower():
                    j += 1
                intron_num += 1
                self.introns.append((intron_num, i, j))
                i = j       # move i to j for new location of next region
            elif seq[i].isupper():
                j = i
                while j < n and seq[j].isupper():
                    j += 1
                exon_num += 1
                self.exons.append((exon_num, i, j))
                i = j
        return (self.introns, self.exons)

# objects
regions = [SplicingRegion(header, seq) for header,seq in records]   # sets up objects for all fasta records
mode = regions[0].mode    
# ^ detects DNA or RNA mode for first object in list of splicingregion objects, 
# only check once b/c all types should be same in the fasta file
# .mode and not .mode() b/c pulling value not running fxn
motif_objects = [Motif(motif, mode) for motif in motifs]


class MotifScanner:
    ''' Scans sequence for motifs and returns identififed locations/hits '''


class MotifLocation:
    ''' Holds record for motif hit 
        Where the motif is found '''
    pass

class MotifMarkRenderer:
    pass