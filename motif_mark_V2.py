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

def color_brew_inator():
    '''Loads the color palette in from a user given file'''
    colors_list = []
    while True:
        line =  COLORS_FILE.readline().rstrip()
        if line == '':
            break
        color_rgb_values = re.findall('[0-9]{1,3}', line)
        colors_list.append(color_rgb_values)
    return colors_list
    #grabbing all of the colors. 
    #   The rgb values of each color are a list. 
    #   All color lists are stored in an outer list (COLORS_LIST)'

COLORS_LIST = color_brew_inator()

def make_colors_ints_inator():
    '''Changes the rgb values in the colors list from strings to ints'''
    colors_ints_list = []
    for color in COLORS_LIST:
        color_ints = []
        for rgb_string in color:
            rgb_int = int(rgb_string)
            color_ints.append(rgb_int)
        colors_ints_list.append(color_ints)
    return colors_ints_list

COLORS_INTS_LIST = make_colors_ints_inator()

print(COLORS_INTS_LIST)

class Gene:
    '''contains header, gene length (total read from fasta, and gene #, method 
    to draw itself'''
    def __init__(self, header, length, gene_count):
        self.header = header
        self.length = length
        self.gene_count = gene_count

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

class Exon:
    '''contains start_pos_pos, gene length (total read from fasta), length,
     and gene # & method to draw itself'''
    def __init__(self, start_pos, length, gene_count):
        self.start_pos = start_pos
        self.length = length
        self.gene_count = gene_count

    def draw(self):
        x = LEFT_MARGIN + self.start_pos
        y = self.gene_count * 30 + 30
        context.set_line_width(20)
        context.move_to(x ,y)
        context.line_to(x+self.length, y)
        context.stroke()

class Motif:
    '''contains start, length, gene_count)'''
    def __init__(self, start_pos, length, gene_count):
        self.start_pos = start_pos
        self.length = length
        self.gene_count = gene_count

    def draw(self):
        x = LEFT_MARGIN + self.start_pos
        y = self.gene_count * 30 + 30
        context.set_line_width(20)
        context.move_to(x ,y)
        context.line_to(x+self.length, y)
        context.stroke()
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
    exon_examp = Exon(100,100,1)
    exon_examp.draw()
beautiful_picture_inator()

