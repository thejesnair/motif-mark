#!/usr/bin/env python

'''
Docstring for motif-mark.motif-mark-oop
Takes in FASTA file and generates image marking motifs along gene
Script outputs both an svg and png
'''

import argparse
import re # regex
from collections import defaultdict
import cairo

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
                if current_header is not None:  # if there is already have a record in store, save
                    records.append((current_header, ''.join(current_seq)))

                # start tracking new record
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)    # reading seq lines for the current record
        
        # Saving last record
        if current_header is not None:
            records.append((current_header, ''.join(current_seq)))
    return records


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
        self.pattern = motif.upper()    # normalize for checking
        if mode == "DNA":
            self.degenerate_map = self.DNA_degenerate    # have to reference class variables
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
        self.regex_body = self.build_regex_body()   # call method
        self.lookahead_overlap = self.build_lookahead_overlap_regex()

    ## METHODS
    def build_regex_body(self) -> str:
        ''' Build regex expression for motif '''
        regex_exp = []
        for char in self.pattern.upper():    # upper bc of dict
            if char in self.degenerate_map:
                regex_exp.append("[" + self.degenerate_map[char] + "]")
            else:
                regex_exp.append(char)
        return ''.join(regex_exp)

    def build_lookahead_overlap_regex(self) -> str:
        ''' Build regex expression for lookahead/overlap '''
        return "(?=(" + self.regex_body + "))"    # (?=...) is a lookahead assertion!

class SplicingRegion:
    ''' Represents FASTA record
        Where I'm looking for motifs '''
    
    def __init__(self, header: str, sequence: str):
        self.header = header
        self.sequence = sequence    # original sequence
        self.sequence_upper = sequence.upper()    # to identify motifs

        self.introns, self.exons = [], []

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
        i = 0    # start coordinate
        j = 0    # 'scanning' end coordinate

        while i < n:
            if seq[i].islower():
                j = i    # start j at current i location
                while j < n and seq[j].islower():
                    j += 1
                intron_num += 1
                self.introns.append((intron_num, i, j))
                i = j    # move i to j for new location of next region
            elif seq[i].isupper():
                j = i
                while j < n and seq[j].isupper():
                    j += 1
                exon_num += 1
                self.exons.append((exon_num, i, j))
                i = j
        return (self.introns, self.exons)


class MotifScanner:
    ''' Scans sequence for motifs and returns identififed locations/hits '''
    
    def __init__(self, region, motifs):
        self.region = region
        self.motifs = motifs

    ## METHODS
    def scan(self):
        ''' Scans sequence for motifs, returns list of MotifLocation '''
        #hits = {}    # dict of key,value pairs motif, hits
        seq = self.region.sequence_upper
        locations = []

        for motif in self.motifs:
            pattern = motif.lookahead_overlap    # pull out regex expression for lookahead assertion
            #hits[motif.pattern] = []    # initialize empty list for each motif pattern

            for match in re.finditer(pattern, seq):    # .finditer(): https://docs.python.org/3/library/re.html#finding-all-adverbs
                start = match.start()
                end = start + motif.length    # comes from self.length in Motif class
                #hits[motif.pattern].append([start, end])    # can append list to empty list value in dict
                location = MotifLocation(
                    self.region.header,    # comes from splicingregion obj passed into motifscanner
                    motif.pattern,    # comes from motif object
                    start,    # from match.start()
                    end    # from start + match.length
                )
                locations.append(location)
        return locations


class MotifLocation:
    ''' Holds record for ONE motif hit
        "Itemize" motif locations, so it's easier to render final image
        Interacts with MotifScanner to hold record
        Where the motif is found '''
    
    # record
    def __init__(self, header, motif_pattern, start, end):
        self.header = header
        self.motif_pattern = motif_pattern
        self.start = start
        self.end = end
        
        self.lane = None    # needed for renderer later on


class MotifMarkRenderer:
    ''' Takes in MotifLocation object and renders image
        Computes lanes for motifs and overlap, renders introns, exons, and motifs '''

    ## METHODS
    def group_headers(self, locations):
        ''' Groups genes/records by header and then stores list of all motif hits
            Sorts hits by start position so that hits are in chronological order left -> right for later lane assignment and drawing '''
        grouped = defaultdict(list)    # initializes an empty list every time a new key is added to dict!

        # first, group headers and their hits
        for loc in locations:
            grouped[loc.header].append(loc)    # stores the full tuple

        # second, sort hits by start position -> useful for lane assignment later on
        for header in grouped:
            grouped[header].sort(key=lambda loc: loc.start)    # (loc.start, loc.end) incase motifs have overlapping starts?
            # lambda: anonymous, small inplace fxn: in this case sorting by start loc, takes loc and returns loc.start, sorts by loc.start
            # some sources:
                # https://www.geeksforgeeks.org/python/python-sort-list-of-list-by-specified-index/
                # https://docs.python.org/3/howto/sorting.html#key-functions
                # https://docs.python.org/3/library/operator.html key=attrgetter("start") use in future?
        return grouped
    
    def assign_lanes(self, hits):
        ''' Takes in MotifLocation list of hits and sorts motif hits by start coordinate
            format: (header, motif, start, end)
            Assigns lane based on coordinates and overlap
            Returns list of lanes and their corresponding hits '''
        lanes = []     # list of lanes: each lane contains a list of hits
        lane_ends = []     # tracking: holds end position for each lane, 0-based

        # loop through all hits, process one at a time
        for hit in hits:
            # assign start and stop
            start = hit.start  # extract coordinates for readability
            end = hit.end

            hit_placed = False  # tracker: is hit assigned to a lane? assume NO
            
            # search for lane to assign hit to
            for lane_pos in range(len(lane_ends)):  # checking every lane that CURRENTLY exists [0, 1, 2...]
                if start >= lane_ends[lane_pos]:    # does hit start AFTER the last hit in this lane ends?
                    lanes[lane_pos].append(hit)     # YES: add hit for that lane position for drawing later
                    lane_ends[lane_pos] = end       # AND update the end coord for that lane
                    hit_placed = True       # update tracker, the hit has been assigned to a lane!
                    break   # stop checking the other lanes, bc hit has been placed into an appropriate lane

            #  if no lane identified: OVERLAP, hit needs to be assigned to a new lane
            if not hit_placed:
                lanes.append([hit])
                lane_ends.append(end)
        return lanes

    # IMAGE RENDERING
    def __init__(self, locations):
        # default for drawing surface
        self.width, self.height = 800, 800
        self.left_margin = 220    # clear space: can place labels in left margin!
        self.right_margin = 40
        self.row_height = 120    # vertical distance b/w records
        self.top_margin = self.row_height    # set equal to row_height bc I want even spacing b/w records
        self.line_width = 3

        # exon
        self.exon_height = 16

        # motifs
        self.motif_height = 15
        self.motif_offset = 25 # above backbone
        self.lane_gap = 14   # vertical space b/w lanes


        # motif color, records
        ''' Note: Only four motifs in assignment so four colors here.
            Future changes to this script will need to potentially account for unknown number of motifs '''
        self.motif_palette = [      # https://rgbcolorpicker.com/0-1
            (0.5, 0.4, 0.6),
            (0.388, 0.533, 0.82),
            (0.96, 0.82, 0.56),
            (0.65, 0.1, 0.25)
        ]

        self.motif_colors = {}

        for loc in locations:
            motif = loc.motif
            if motif not in self.motif_colors:
                index = len(self.motif_colors) #% len(self.motif_palette)    # if uneven num of motifs + colors, assign pos
                self.motif_colors[motif] = self.motif_palette[index]


        # header labels
        self.font_size = 16
        self.label_x = 80    # position for label in left margin

        # create PNG surface
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)
        self.ctx = cairo.Context(self.surface)

        # background
        self.ctx.set_source_rgb(1,1,1)    # need to set white background, default is transparent
        self.ctx.paint()

        # LEGEND #
        self.legend_x = self.width - self.right_margin - 80
        self.legend_y = 30
        self.legend_box_size = 14
        self.legend_gap = 22

        # font
        self.ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)    # set font type
        self.ctx.set_font_size(self.font_size)

        # intron key
        self.line_y = self.legend_y + (self.legend_box_size/2)
        self.ctx.set_line_width(self.line_width)
        self.ctx.move_to(self.legend_x, self.line_y)
        self.ctx.line_to(self.legend_x + self.legend_box_size, self.line_y)
        # draw line
        self.ctx.set_source_rgb(0,0,0)    # need to set line color to black
        self.ctx.set_line_width(3)
        self.ctx.stroke()

        self.ctx.move_to(self.legend_x + self.legend_box_size + 10, self.legend_y + self.legend_box_size)
        self.ctx.show_text("Intron")

        self.legend_y += self.legend_gap

        # exon
        self.ctx.rectangle(self.legend_x, self.legend_y, self.legend_box_size, self.legend_box_size)
        self.ctx.fill()

        self.ctx.move_to(self.legend_x + self.legend_box_size + 10, self.legend_y + self.legend_box_size)
        self.ctx.show_text("Exon")

        self.legend_y += self.legend_gap


    def render_motifs(self, region, locations):

        # draw legend box
        for motif_name, (r,g,b) in motif_colors.items():
            ctx.set_source_rgb(r,g,b)
            ctx.rectangle(legend_x, legend_y + 20, legend_box_size, legend_box_size)
            ctx.fill()

            # text labels
            ctx.set_source_rgb(0,0,0)
            ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) # this will apply to headers if not changed
            ctx.set_font_size(font_size)    # will also apply to headers if not changed
            ctx.move_to(legend_x + legend_box_size +10, legend_y+20 + legend_box_size)
            ctx.show_text(motif_name)

            legend_y += legend_gap + 20

        for i, (header, length_bp, exons, motifs) in enumerate(records):
            y = top_margin + (i * row_height)    # space between each record vertically (lane stacking)
            # i: for each record * row_gap will create even gaps

            x0 = left_margin    # record starting pos
            x1 = left_margin + length_bp    # genomic length/final ending position

            # drawn backbone
            ctx.move_to(x0, y)
            ctx.line_to(x1, y)
            ctx.stroke()

            # header label
            ctx.set_source_rgb(0,0,0)
            #ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            #ctx.set_font_size(font_size)
            
            ctx.move_to(label_x, y+5)
            ctx.show_text(header)

            # exons
            for (exon_start, exon_end) in exons:
                exon_x = left_margin + exon_start
                exon_w = exon_end - exon_start  # how long is exon
                exon_y = y - (exon_height/2)

                ctx.rectangle(exon_x, exon_y, exon_w, exon_height)  # (x0, y0, x1, y1)
                ctx.set_source_rgb(0,0,0)   # black fill for now
                ctx.fill()

            # motifs
            for (lane_pos, lane_hits) in enumerate(motifs):
                lane_y = y - motif_offset - lane_pos * lane_gap
            
                for (motif_name, motif_start, motif_end) in lane_hits:
                    motif_x = left_margin + motif_start
                    motif_w = motif_end - motif_start
                    motif_y = y - motif_offset

                    # draw motif in specified color
                    r,g,b = motif_colors[motif_name]
                    ctx.set_source_rgb(r,g,b)   # draw motif
                    ctx.rectangle(motif_x, motif_y, motif_w, motif_height)
                    ctx.fill()
                    ctx.set_source_rgb(0,0,0)   # change back to black for next record

        surface.write_to_png("motif_mark_test4.png")










def main():

    def get_args():
        '''Takes in command line arguments for paths to motifs files, FASTA file, and output image.'''
        parser = argparse.ArgumentParser(description=
                                        "Program to generate an image of motifs along a gene")
        parser.add_argument("-f", "--file", help="path for fasta file", type=str, required=True)
        parser.add_argument("-m", "--motif", help="path for file of all possible motifs", type=str, required=True)
        #parser.add_argument("-o", "--outfile", help = "name for outputted motif image (png)", type=str, required=True)

        return parser.parse_args()

    args = get_args()
    input_fasta = args.file
    motif_file = args.motif
    #image = args.outfile
    
    records = read_fasta(input_fasta)

    # motifs list
    with open(motif_file) as f:
        motifs = [line.strip() for line in f]   # strips trailing space, adds motifs to list

    print(motifs)   #check that properly pulls motifs

    # objects needed for MotifScanner
    regions = [SplicingRegion(header, seq) for header,seq in records]   # sets up objects for all fasta records
    mode = regions[0].mode    
    # ^ detects DNA or RNA mode for first object in list of splicingregion objects, 
    # only check once b/c all types should be same in the fasta file
    # .mode and not .mode() b/c pulling value not running fxn
    motif_objects = [Motif(motif, mode) for motif in motifs]

    # region object contains: uppercase sequence, sequence, header
    # motif object contains: lookahead regex, length, regex


if __name__ == "__main__":
    #main()

    # TESTING
    test_region = SplicingRegion(">test", "aaGGcccTTTaGGGGG")
    print(test_region.header)
    print(test_region.sequence_upper)
    print(test_region.mode)

    mode = test_region.mode
    m1 = Motif("YGCY", mode)

    print(m1.pattern)
    print(m1.lookahead_overlap)
    print(m1.length)
    scanner = MotifScanner(test_region, [m1])

    print(scanner.region.header)
    print(len(scanner.motifs))

    print(scanner.scan())

    render = MotifMarkRenderer()