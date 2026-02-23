#!/usr/bin/env python

"""
testing.py is a 'scratchwork' script used to do standalone testing of smaller functions and
write quick pseudocode as I develop classes/functions

"""

# pattern = "YYYYYYYYYY"


# DNA_degenerate = {'Y': "CT", 'R': "AG", 'N': "ACGT",
#                     'V': "ACG", 'H': 'ACT', 'D': "AGT",
#                     'B': "CGT"}

# ## METHODS
# def build_regex_body(pattern, degenerate_map):
#     regex_exp = []
#     for char in pattern.upper():
#         if char in degenerate_map:
#             regex_exp.append("[" + degenerate_map[char] + "]")
#         else:
#             regex_exp.append(char)
#     return ''.join(regex_exp)


# print(build_regex_body(pattern, DNA_degenerate))


"""
Note for later:
keep regex upper case, when scanning sequence--scan uppercase but keep original
then can map motif matches onto original sequence with position, that way wont miss motifs in either intron or exon regions"""


# def read_fasta(input_fa:str) -> str:
#     ''' Parses a FASTA file and returns a list of (header, sequence) tuples
#         Handles multiline sequences and stores seq as oneline '''
    
#     current_header = None
#     current_seq = []
#     records = []
    
#     with open(input_fa, 'r') as in_fa:
#         for i, line in enumerate(in_fa):
#             line = line.strip()
#             if line.startswith('>'):    # new header, new record starting
#                 if current_header is not None:      # if there is already have a record in store, save
#                     records.append((current_header, ''.join(current_seq)))

#                 # start tracking new record
#                 current_header = line
#                 current_seq = []
#             else:
#                 current_seq.append(line)        # reading seq lines for the current record
        
#         # Saving last record
#         if current_header is not None:
#             records.append((current_header, ''.join(current_seq)))
#     return records

# print(read_fasta("/Users/thejesnair/Desktop/KCGIP/Bi625/Motif_mark/motif-mark/Figure_1.fasta"))


# seq = "atgtccacatgtagtcacgtttgacatcccagggccacctcagcaggccgtctctgggga"
# seq2 = "AUGUCCACAUGUAGUCAGCUUU"
# seq3 = "ilovecake"

# def detect_mode(seq) -> str:
#     ''' Determine if DNA or RNA '''
#     if "U" in seq.upper():
#         return "RNA"
#     elif "T" in seq.upper():
#         return "DNA"
#     else:
#         raise ValueError("DNA or RNA FASTA required")
    

#print(detect_mode(seq))
#print(detect_mode(seq2))
#print(detect_mode(seq3))


# seq = "aaBBcccDDDDeFFFFF"
# introns = []
# exons = []
# def find_intron_exon_regions(seq):
#     ''' Determines coordinates of introns and exons
#         Tuple: Intron/exon number, start pos, end pos '''

#     intron_num = 0
#     exon_num = 0
#     n = len(seq)
#     i = 0       # start coordinate
#     j = 0       # 'scanning' end coordinate

#     while i < n:
#         if seq[i].islower():
#             j = i       # start j at current i location
#             while j < n and seq[j].islower():
#                 j += 1
#             intron_num += 1
#             introns.append((intron_num, i, j))
#             i = j       # move i to j for new location of next region
#             j = 0       # reset j so loop can continue sliding
#         elif seq[i].isupper():
#             j = i
#             while j < n and seq[j].isupper():
#                 j += 1
#             exon_num += 1
#             exons.append((exon_num, i, j))
#             i = j
#             j = 0
#     return (introns, exons)


# print(find_intron_exon_regions(seq))


'''
variables:
list of intron coordinates
list of exon coordinates

function
variables:
intron num
exon num
start coordinate
end coordinate
length of input sequence

loop logic

while the start position < length of sequence:
if the char in sequence is lower and the start coordinate is < length of sequence
set end coordinate to plus one of start coordinate
while ... and lowercase, keep adding 1 to end coordiante
need to reset end coordinate and shfit start coordinate for every new intron/exon encountered
if reaches uppercase then set j coordiante to new position and continue loop
return intron/exon list of tuples
'''


'''

class MotifScanner

def __init__(self, splicing_object, motif_object):
    sequence
    sequence.upper
    regex
    regex.lookahead

    def fxn
    scan sequence for regex(upper), scan sequence for regex.lookahead

    return list of all motifs and their hits (position/indexing)
    
    '''


# import re

# seq = "AATGCCTGCTAGGGGG"
# motif = "YGCY"
# pattern = "(?=([CT]GC[CT]))"

# hits = {}
# for m in re.finditer(pattern, seq):
#     start = m.start()
#     motif_len = len(motif)
#     end = start + motif_len
#     hits[motif] = [start, end]

# print(hits)


#Motif mark renderer

'''
assigning motifs to lanes if there's overlap
lanes arent necessarily per motif, they're for overlapping PER INTERVAL

so for one gene:
introns are just baseline
exons a dark grey/black box
each motif assigned a color
assign lanes based on overlap, otherwise if a single motif is within an interval theyre drawn at a default position
lane 0 (default)
lane 1: if overlaps with 0
lane 2: if overlap with 1 and 0 

'''

'''
1. group motif hits by gene/header
2. create cairo image for ONE output image
3. for each fasta record create a baseline drawing (stack them vertically)
    a. get hits for each record
    b. assign lanes (and overlaps) for record's hits
    c. draw gene "backbone": introns, exons
    d. draw motifs using assigned lanes
4. draw motif KEY (outside of loop!)
5. finish and output image
'''

loc1 = ("header1", "YGCY", 2, 6)
loc2 = ("header1", "CATAG", 10, 15)
loc3 = ("header1", "YGCY", 4, 8)
locations = [loc1, loc2, loc3]

# def group_headers(locations):
#     grouped = {}
#     for loc in locations:
#         header = loc[0]
#         if header not in grouped:
#             grouped[header] = []
#         grouped[header].append(loc)     # stores the full tuple
#     return grouped
        
# print(group_headers(locations))


'''
def assign_lanes()
take in locations (dont have to worry about sorting by header/gene because each object will be assigned to one fasta record)
sort locations by (start,end) [shouldnt it already be sorted?]
    - no theyre not already sorted bc it depends on the motif!

    look at start and end coordinates to determine overlap
    ex: (motifhit1, 2, 45) vs (motifhit2, 10, 55): overlap from (10,45)

    lane_ends = [end position of latest hit per lane]

    exs:
    hit1 (2,45)
    hit2 (10,55)
    hit3 (35,65)
    hit4 (55, 80)

    1. sort hits by start coordinate
        if a tie then sort by end coordinate
    2. hit1 -> no hits exist so default to lane 0
       update lane_ends = [45]
    3. hit2 -> 10 >=45, false
       make lane 1
       update lane 1 end to 55
       lane_ends = [45,55]
    4. hit3 -> 35 >=45, false
       35 >= 55, false
       update lane 2 end to 65
       lane_ends = [45,55,65]
    5. hit4 -> 55 >= 45, true
       add to lane 0
       update lane 0 end to 55
       lane_ends = [80, 55, 65]
'''

'''
for each hit:
    get start/end coordinate
    assume hit is NOT placed yet
    
    for each lane that EXISTS:
        if hit fits in lane:
                add hit to the lane
                update lane end coord
                mark as assigned
                and then stop checking lanes?
        if hit has never been assigned:
            create a new lane
            assign hit
            update end coord
            
'''

# def assign_lanes(locations):
#     # locations: list of hits sorted by start coord
#     # hit format: (header, motif, start, end)
#     lanes = []     # list of lanes: each lane contains a list of hits
#     lane_ends = []     # tracking: holds end position for each lane, 0-based

#     # loop through all hits, process one at a time
#     for hit in locations:
#         # assign start and stop
#         start = hit[2]  # extract coordinates for readability
#         end = hit[3]

#         hit_placed = False  # tracker: is hit assigned to a lane? assume NO
        
#         # search for lane to assign hit to
#         for lane_pos in range(len(lane_ends)):  # checking every lane that CURRENTLY exists [0, 1, 2...]
#             if start >= lane_ends[lane_pos]:    # does hit start AFTER the last hit in this lane ends?
#                 lanes[lane_pos].append(hit)     # YES: add hit for that lane position for drawing later
#                 lane_ends[lane_pos] = end       # AND update the end coord for that lane
#                 hit_placed = True       # update tracker, the hit has been assigned to a lane!
#                 break   # stop checking the other lanes, bc hit has been placed into an appropriate lane

#         #  if no lane identified: OVERLAP, hit needs to be assigned to a new lane
#         if not hit_placed:
#             lanes.append([hit])
#             lane_ends.append(end)
#     return lanes

    
# print(assign_lanes(locations))


'''

PYCAIRO IMAGE RENDERER
1. Draw backbone (length of record)
2. Drawn exons as boxes on line (non boxed lines will represent introns)
3. motifs + lanes
4. need to draw key, label each drawing with header
5. output image as PNG for submission, not SVG

'''
import cairo

# default settings for surfacce
width, height = 800, 800
left_margin = 220    # clear space: can place labels in left margin!
right_margin = 40
row_height = 120    # vertical distance b/w records
top_margin = row_height    # set equal to row_height bc I want even spacing b/w records
line_width = 3

# exon
exon_height = 16

# motifs
motif_height = 15
motif_offset = 25 # above backbone
lane_gap = 14   # vertical space b/w lanes

# motif color, records
motif_colors = {
    "YGCY": (.5,.4,.6),
    "CATAG": (0.388, 0.533, 0.82),
    "YYYY": (0.96, 0.82, 0.56)
}

records = [
    ("header1", 200, 
     [(20,60), (120,160)],  # exon
     [                      # motifs by lane
        [("YGCY", 30,40), ("CATAG", 130, 145)],  # lane 0
        [("YGCY", 35, 55)]              # lane 1 (overlaps lane 0)
        ]),
    ("header2", 400, 
     [(50,100), (200,260), (330,380)],   # exons
     [
     [("YYYY", 45, 68), ("CATAG", 112, 160)]    # lane 0
        ]),
    ("header3", 340, 
     [(32, 70), (160, 210)], # exons
     [
     [("YGCY", 78, 98), ("CATAG", 132, 145)]    # lane 0
        ])
]

# header labels
font_size = 16
label_x = 80    # position for label in left margin

# create PNG surface
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
ctx = cairo.Context(surface)

# background
ctx.set_source_rgb(1,1,1)    # need to set white background, default is transparent
ctx.paint()

# draw line
#ctx.set_source_rgb(0,0,0)    # need to set line color to black
#ctx.set_line_width(3)
#ctx.move_to(100,58)    # (x,y) (0,0 is top left) x: right, y: down
#ctx.line_to(420,58)
#ctx.stroke()

drawing_margin = width - left_margin - right_margin

# legend
legend_x = width - right_margin - 80
legend_y = 30
legend_box_size = 14
legend_gap = 22
#legend_font_size = 14

ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
ctx.set_font_size(font_size)


# intron key
line_y = legend_y + (legend_box_size/2)
ctx.set_line_width(line_width)
ctx.move_to(legend_x, line_y)
ctx.line_to(legend_x + legend_box_size, line_y)
# draw line
ctx.set_source_rgb(0,0,0)
ctx.set_line_width(3)
ctx.stroke()

ctx.move_to(legend_x + legend_box_size + 10, legend_y + legend_box_size)
ctx.show_text("Intron")

legend_y += legend_gap

# exon
ctx.rectangle(legend_x, legend_y, legend_box_size, legend_box_size)
ctx.fill()

ctx.move_to(legend_x + legend_box_size + 10, legend_y + legend_box_size)
ctx.show_text("Exon")

legend_y += legend_gap

# draw legend box
for motif_name, (r,g,b) in motif_colors.items():
    ctx.set_source_rgb(r,g,b)
    ctx.rectangle(legend_x, legend_y, legend_box_size, legend_box_size)
    ctx.fill()

    # text labels
    ctx.set_source_rgb(0,0,0)
    ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) # this will apply to headers if not changed
    ctx.set_font_size(font_size)    # will also apply to headers if not changed
    ctx.move_to(legend_x + legend_box_size +10, legend_y + legend_box_size)
    ctx.show_text(motif_name)

    legend_y += legend_gap

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
            motif_y = lane_y

            # draw motif in specified color
            r,g,b = motif_colors[motif_name]
            ctx.set_source_rgb(r,g,b)   # draw motif
            ctx.rectangle(motif_x, motif_y, motif_w, motif_height)
            ctx.fill()
            ctx.set_source_rgb(0,0,0)   # change back to black for next record

surface.write_to_png("motif_mark_test5.png")

