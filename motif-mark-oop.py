#!/usr/bin/env python

'''
motif-mark.py generates an image marking motifs along a gene utilizing object oriented programming
Script takes in a FASTA file and txt file of motifs and outputs a PNG

Requirements: pycairo installation (conda install -n my_pycairo pycairo)
'''

import argparse
import re   # regex
from collections import defaultdict
import cairo
from pathlib import Path    # so I can extract name for output image

## ARGPARSE ##
def get_args():
    ''' Takes in command line arguments for paths to motifs file and FASTA file
        Note: Script requires pycairo installation to run '''
    parser = argparse.ArgumentParser(description=
                                    "Program to generate a PNG of motifs along a gene")
    parser.add_argument("-f", "--file", help="path for FASTA file", type=str, required=True)
    parser.add_argument("-m", "--motif", help="path for txt file of possible motifs", type=str, required=True)

    return parser.parse_args()

## HELPER FXN ##
# read in fasta
def read_fasta(input_fa:str) -> list[tuple[str,str]]:
    ''' Parses FASTA file, handles multiline sequences and stores seq as oneline
        Returns: 
            list[(header, sequence)] '''
    
    current_header = None
    current_seq = []
    records = []
    
    with open(input_fa, 'r') as in_fa:
        for line in in_fa:
            line = line.strip()
            if line.startswith('>'):    # new header, new record starting
                if current_header is not None:  # if there is already a record in store, save
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
        # https://stackoverflow.com/questions/11430863/how-to-find-overlapping-matches-with-a-regexp

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
            Returns for introns and exons:
                list[(number, start pos, end pos)] '''

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
        ''' Scans sequence for motifs
            Uses MotifLocation class to hold record
            Returns:
                list[MotifLocation(header, motif pattern, start coord, end coord)] '''
        seq = self.region.sequence_upper
        locations = []

        for motif in self.motifs:
            pattern = motif.lookahead_overlap    # pull out regex expression for lookahead assertion

            for match in re.finditer(pattern, seq):    # .finditer(): https://docs.python.org/3/library/re.html#finding-all-adverbs
                start = match.start()
                end = start + motif.length    # comes from self.length in Motif class
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
        Computes lanes for motifs and overlap, renders introns, exons, and motifs
        Returns:
            Hits, sorted by header and chronological order
            Lane assignments for hits '''

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
            Returns:
                List of lanes and corresponding hits '''
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
    def __init__(self, regions, locations):
        # default for drawing surface
        n_records = len(regions)
        self.width = 1200
        self.left_margin = 220    # clear space: can place labels in left margin!
        self.right_margin = 40
        self.row_height = 200    # vertical distance b/w records
        self.top_margin = self.row_height    # set equal to row_height bc I want even spacing b/w records
        self.line_width = 3
        
        # exon
        self.exon_height = 16

        # motifs
        self.motif_height = 15
        self.motif_offset = 25 # above backbone
        self.lane_gap = 14   # vertical space b/w lanes

        # bottom margin + height computed from last record, implemented so image can render more records
        self.bottom_margin = 40

        # y position of last record
        if n_records == 0:
            last_y = self.top_margin
        else:
            last_y = self.top_margin + (n_records - 1) * self.row_height

        # canvas height
        self.height = int(last_y + (self.exon_height/2) + self.bottom_margin)

        # motif color, records
        ''' Note: Up to five motifs supported '''
        self.motif_palette = [      # https://rgbcolorpicker.com/0-1
            (0.5, 0.4, 0.6),        # purple
            (0.30, 0.60, 0.85),     # blue
            (0.87, 0.68, 0.33),     # orangey
            (0.36, 0.63, 0.58),     # teal
            (0.78, 0.52, 0.60)      # rose
        ]

        self.motif_colors = {}      # motifs and assigned colors

        for loc in locations:
            motif = loc.motif_pattern
            if motif not in self.motif_colors:      # if color NOT assigned
                index = len(self.motif_colors) #% len(self.motif_palette)    # commented out code: if uneven num of motifs + colors, assign pos
                self.motif_colors[motif] = self.motif_palette[index]


        # header labels
        self.font_size = 16
        self.label_x = 80    # position for label in left margin

        # create PNG surface
        self.scale = 2  # scale up because image seemed small since motifs are small
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width*self.scale, self.height*self.scale)   # for png output
        self.ctx = cairo.Context(self.surface)
        self.ctx.scale(self.scale, self.scale)

        # background
        self.ctx.set_source_rgb(1,1,1)    # need to set white background, default is transparent
        self.ctx.paint()

        # LEGEND
        self.legend_x = self.width - self.right_margin - 200
        self.legend_y = 30
        self.legend_box_size = 14
        self.legend_gap = 22

        # font
        self.ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)    # set font type
        self.ctx.set_font_size(self.font_size)

        legend_y = self.legend_y    # keep as starting position, will need to rewrite after every loop

        # intron key
        line_y = legend_y + (self.legend_box_size/2)
        self.ctx.set_line_width(self.line_width)
        self.ctx.move_to(self.legend_x, line_y)
        self.ctx.line_to(self.legend_x + self.legend_box_size, line_y)
        # draw line
        self.ctx.set_source_rgb(0,0,0)    # need to set line color to black
        self.ctx.set_line_width(3)
        self.ctx.stroke()

        self.ctx.move_to(self.legend_x + self.legend_box_size + 10, legend_y + self.legend_box_size)
        self.ctx.show_text("Intron")

        legend_y += self.legend_gap

        # exon
        self.ctx.rectangle(self.legend_x, legend_y+4, self.legend_box_size, self.legend_box_size)
        self.ctx.fill()
        # exon key
        self.ctx.move_to(self.legend_x + self.legend_box_size + 10, legend_y + self.legend_box_size + 4)
        self.ctx.show_text("Exon")

    def render_motifs(self, regions, locations, output_png):
        grouped = self.group_headers(locations)     # header contains list[MotifLocation]
        
        legend_y = self.legend_y  + 2.4 * self.legend_gap   # space b/w exon label and motif labels
        
        # draw legend box
        for motif_name, (r,g,b) in self.motif_colors.items():
            self.ctx.set_source_rgb(r,g,b)
            self.ctx.rectangle(self.legend_x, legend_y, self.legend_box_size, self.legend_box_size)
            self.ctx.fill()

            # text labels
            self.ctx.set_source_rgb(0,0,0)
            self.ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD) # this will apply to headers if not changed
            self.ctx.set_font_size(self.font_size)    # will also apply to headers if not changed
            self.ctx.move_to(self.legend_x + self.legend_box_size +10, legend_y + self.legend_box_size)
            self.ctx.show_text(motif_name)

            legend_y += self.legend_gap + 5     # space between motif labels

        for i, region in enumerate(regions):
            y = self.top_margin + (i * self.row_height)    # space between each record vertically (lane stacking)
            # i: for each record * row_height will create even gaps

            length_bp = len(region.sequence)
            x0 = self.left_margin    # record starting pos
            x1 = self.left_margin + length_bp  # genomic length/final ending position

            # draw backbone
            self.ctx.move_to(x0, y)
            self.ctx.line_to(x1, y)
            self.ctx.stroke()

            # header label
            self.ctx.set_source_rgb(0,0,0)

            self.ctx.move_to(self.label_x, y-140)   # HEADER SPACING FROM RECORD
            self.ctx.show_text(region.header)

            # exons
            for (_, exon_start, exon_end) in region.exons:  # _ means intentionally not using this value! dont need exon_num for rendering
                exon_x = self.left_margin + exon_start
                exon_w = exon_end - exon_start  # how long is exon?
                exon_y = y - (self.exon_height/2)

                self.ctx.rectangle(exon_x, exon_y, exon_w, self.exon_height)  # (x0, y0, x1, y1)
                self.ctx.set_source_rgb(0,0,0)   # black fill for now
                self.ctx.fill()

            # motifs -> lanes
            hits = grouped.get(region.header, [])
            lanes = self.assign_lanes(hits)

            for (lane_pos, lane_hits) in enumerate(lanes):
                lane_y = y - self.motif_offset - lane_pos * self.lane_gap
            
                for hit in lane_hits:
                    motif_name = hit.motif_pattern
                    motif_x = self.left_margin + hit.start
                    motif_w = hit.end - hit.start
                    motif_y = lane_y    # bug in dev script: LANE ASSIGNMENT!!! without this other lanes aren't drawn

                    # draw motif in specified color
                    r,g,b = self.motif_colors[motif_name]
                    self.ctx.set_source_rgb(r,g,b)   # draw motif
                    self.ctx.rectangle(motif_x, motif_y, motif_w, self.motif_height)
                    self.ctx.fill()
                    self.ctx.set_source_rgb(0,0,0)   # change back to black for next record

        self.surface.write_to_png(output_png)


def main():

    # argparse/file vars
    args = get_args()   
    input_fasta = args.file
    motif_file = args.motif
    output_png = Path(input_fasta).with_suffix(".png") # Path pulls exact name from input fasta for output png

    records = read_fasta(input_fasta)

    # motifs list
    with open(motif_file) as f:
        motif_strings = [line.strip() for line in f]   # strips trailing space, adds motifs to list

    # objects needed for MotifScanner
    regions = [SplicingRegion(header, seq) for header,seq in records]   # sets up objects for all fasta records
    mode = regions[0].mode  # identify mode for all records, just need to check first record

    # normalize motifs (needs to happen AFTER mode assignment)
    if mode == "DNA":
        motif_strings = [m.upper().replace("U", "T") for m in motif_strings]
    elif mode == "RNA":
        motif_strings = [m.upper().replace("T", "U") for m in motif_strings]

    motif_objects = [Motif(motif, mode) for motif in motif_strings]

    # region object contains: uppercase sequence, sequence, header
    # motif object contains: lookahead regex, length, regex

    # scan everything! and put it all in MotifLocation hits
    all_locations = []
    for region in regions:
        scanner = MotifScanner(region, motif_objects)
        hits = scanner.scan()   # return list of MotifLocation
        all_locations.extend(hits)  # adds elements individually (extend flattens), append would nest them

    renderer = MotifMarkRenderer(regions, all_locations)
    renderer.render_motifs(regions, all_locations, output_png)

if __name__ == "__main__":
    main()
