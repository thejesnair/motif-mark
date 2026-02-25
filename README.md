# Motif Mark OOP

motif-mark.py utilizes object-oriented programming to generate a PNG of motifs along gene regions from a FASTA file. Exons, introns, and motif hits are rendered using Pycairo and are to scale. Motifs are color coded and overlapping hits are staggered through automatic lane assignment. Python script currently is able to handle RNA or DNA records, degenerate bases, corrects for base errors in motifs (ex: U instead of T in DNA records). This script can visualize up to 4 motifs.

In order to run, **motif-mark.py requires installation of a pycairo environment**. 

The outputted png retains the same file name as the inputted FASTA file.

Ex:<br>
Input: Figure_1.fasta<br>
Output: Figure_1.png

## Classes utilized

- Motif
    - Motif information: motif definition, regex pattern, ambiguity rules
- SplicingRegion
    - FASTA record: DNA/RNA mode, exon/intron coordinates
- MotifScanner
    - Scans for motifs and records coordinates
- MotifLocation
    - Holds information for a single motif hit
- MotifMarkRenderer
    - Renders final image utilizing objects from SplicingRegion and MotifLocation

----------
<br>

## Future development
### Planned improvements:
1. Motif hit summary report
2. Clean header for output image (gene name instead of whole fasta header)
3. SVG visualization (as an alternative or optional output alongside PNG)

### Possible future considerations
1. Alignment gaps
2. Support for 4+ motifs with automatic color generation (human readable colors)
3. Support for additional/longer fasta records with automatic pycairo surface resizing