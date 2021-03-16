#!/usr/bin/python

# Motif Marker V2
# Tucker Bower
# 3/11/21

# This program searches a .fasta file for motifs of interest and produces a vector image graphic illustrating the locations of each given motif
# and exons for each read or gene of that .fasta file to scale.

# Bash run command:
# python motif_mark_V2.py -f example_data/Figure_1.fasta -m example_data/Fig_1_motifs.txt

import argparse
import re
import cairo

parser = argparse.ArgumentParser(description = 'Input: fasta file with exons capitalized and motifs file | Output: Motif location visualization ')
parser.add_argument('-f' , '--fasta' , type = str, nargs = 1, help = 'input .fasta file and file path (ex: ./example_data/Figure_1.fasta')
parser.add_argument('-m' , '--motifs' , type = str, nargs = 1, help = 'file containing motifs of interest, one per line')
parser.add_argument('-o' , '--output' , type = str, nargs = 1, help = 'Output directory', default = './example_data')
parser.add_argument('-c' , '--colors' , type = str, nargs = 1, help = 'Colors file', default = ['./example_data/pastels.txt'])
args = parser.parse_args()

FASTA_FILE = open(args.fasta[0], "r")
MOTIFS_FILE = open(args.motifs[0], "r")
COLORS_FILE = open(args.colors[0], "r")

OUTPUT_DIR = args.output
FASTA_PREFIX = re.findall('/(.+)\..+', args.fasta[0])[0]

OUTPUT_FILE_NAME = OUTPUT_DIR + '/' + FASTA_PREFIX + '.svg'

IUPAC_CODES_DICT = {'A':'[Aa]', 'T':'[TtUu]', 'U':'[TtUu]', 'C':'[Cc]', 'G':'[Gg]', 'C':'[Cc]', 'R':'[AaGg]', 'Y':'[CcTtUU]',
    'W':'[AaTtUu]', 'S':'[GgCc]', 'M':'[AaCc]', 'K':'[GgTtUu]', 'B':'[CcGgTtUu]', 'D':'[AaGgTtUu]', 'H':'[AaCcTtUu]',
    'V':'[AaCcGg]', 'N':'[GgCcAaTtUu]'}
#Dictionary where the keys are IUPAC codes and values are sequence bases

### global variables for pycairo ###
LEFT_MARGIN = 10


class Gene:
    
    def __init__(self, header, length, gene_count):
        '''1 gene = 1 'read' from the .Fasta file, contains introns, exons, motifs'''

        ### Data ###

        self.header = header
        self.length = length
        self.gene_count = gene_count

        ### Methods ### 

    def draw(self):
        x = LEFT_MARGIN
        y = self.gene_count * 30
        context.move_to(x,y)
        context.show_text(self.header)
        y += 30
        context.move_to(x,y)
        context.set_line_width(1)
        context.line_to(x + self.length, y)

        context.stroke()


        #context.show_text(gene)

#
### PLACEHOLDER CODE ###
#
max_seq_length = 800
number_of_genes = 4
#
### PLACEHOLDER CODE ###
#

### initialize pycairo 'context' object (the surface we draw on)
figure_surface = cairo.SVGSurface(OUTPUT_FILE_NAME, max_seq_length + 100, number_of_genes * 75 + 100)
context = cairo.Context(figure_surface)
context.set_source_rgba(0, 0, 0, 1)
context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL)
context.set_font_size(12)

def beautiful_picture_inator():
    gene_examp = Gene('foo', 800, 1)
    gene_examp.draw()

beautiful_picture_inator()

