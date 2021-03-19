#!/usr/bin/python

# Motif Marker V2
# Tucker Bower
# 3/11/21

# This program searches a .fasta file for motifs of interest and produces a vector image graphic illustrating the locations of each given motif
# and exons for each read or gene of that .fasta file to scale.

# Bash run command:
# python motif_mark_V2.py -f example_data/Figure_1.fasta -m example_data/Fig_1_motifs.txt

import argparse
import regex as re
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
TOP_MARGIN = 10
SPACE_BETWEEN_GENES = 60
SPACE_BTW_HEADER_AND_GENE = 30
SPACE_BTW_LEGEND_MOTIFS = 20

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

def fasta_string_inator(fasta_file):
    '''Reads fasta file to create a dictionary where
    keys: gene names
    values: DNA/RNA sequence as string'''
    gene_sequence_dict = {}
    while True:
        line = fasta_file.readline().rstrip()
        if line == '':
            break
            #End of file
        if line.startswith('>'):
            #Header
            header = re.findall('>(.+)' , line)[0]
            #Extracting gene ID
        else:
            if header not in gene_sequence_dict:
                gene_sequence_dict[header] = line
                #Initializing string in dictionary
            else:
                gene_sequence_dict[header] += line
                #Adding all lines of sequence together into a single continuous string 
    return gene_sequence_dict

GENE_SEQUENCE_DICT = fasta_string_inator(FASTA_FILE)

def iupac_regex_inator(motifs_file):
    '''Reads motifs file and creates a dictionary where:
    keys: motifs 
    values: regex terms that can be used to locate that motif, accounting for IUPAC ambiguous nucleotide notation'''
    motif_regex_dict = {}
    while True:
        motif = motifs_file.readline().rstrip()
        motif_regex = ''
        if motif == '':
            break
        #EOF
        for character in motif:
            motif_regex += IUPAC_CODES_DICT[character.upper()]
        #Translating every single NUC character into a regex term for every possible character that COULD be that character
        motif_regex_dict[motif] = motif_regex
        #Storing those regex terms in a dictionary
    return motif_regex_dict

MOTIF_REGEX_DICT = iupac_regex_inator(MOTIFS_FILE)

def motif_coordinates_inator(motif, sequence):
    '''Uses regular expressions to search for motifs in a sequence of bases and returns a list of startpoints and a list of lengths of each occurence of that motif in a list.'''
    motif_starts_list = []
    motif_lengths_list = []
    for motif_match in re.finditer(motif, sequence, overlapped=True):
        motif_start = motif_match.start() + 1
        motif_end = motif_match.end()
        motif_length =  motif_end - motif_start +1
        motif_starts_list.append(motif_start)
        motif_lengths_list.append(motif_length)
    return motif_starts_list, motif_lengths_list

# test, west = motif_coordinates_inator ('tag','atagatag')
# print('starts', test)
# print(' ')
# print('lengths', west)

class Header:
    '''contains the header, draws the header'''
    def __init__(self, header, gene_number):
        self.header = header
        self.gene_number = gene_number

    def draw(self):
        x = LEFT_MARGIN
        y = self.gene_number * SPACE_BETWEEN_GENES + TOP_MARGIN
        context.move_to(x,y)
        context.set_source_rgba(0, 0, 0, 1)
        context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        context.set_font_size(12)
        context.show_text(self.header)
        context.stroke()

class Gene:
    ''' gene length (total read from fasta, and gene #, method 
    to draw itself'''
    def __init__(self, length, sequence, gene_number):
        self.length = length
        self.sequence = sequence
        self.gene_number = gene_number

    def draw(self):
        x = LEFT_MARGIN
        y = self.gene_number * SPACE_BETWEEN_GENES + SPACE_BTW_HEADER_AND_GENE + TOP_MARGIN
        context.move_to(x,y)
        context.set_source_rgba(0, 0, 0, 1)
        context.set_line_width(1)
        context.line_to(x + self.length, y)
        context.stroke()
        #context.show_text(gene)

class Exon:
    '''contains start_pos_pos, gene length (total read from fasta), length,
     and gene # & method to draw itself'''
    def __init__(self, start_pos, length, gene_number):
        self.start_pos = start_pos
        self.length = length
        self.gene_number = gene_number

    def draw(self):
        x = LEFT_MARGIN + self.start_pos
        y = self.gene_number * SPACE_BETWEEN_GENES + SPACE_BTW_HEADER_AND_GENE +TOP_MARGIN
        context.set_line_width(20)
        context.set_source_rgba(0, 0, 0, 1)
        context.move_to(x ,y)
        context.line_to(x+self.length, y)
        context.stroke()

class Motif:
    '''contains start, length, gene_number)'''
    def __init__(self, start_pos, length, color, gene_number):
        self.start_pos = start_pos
        self.length = length
        self.gene_number = gene_number
        self.color = color

    def draw(self):
        x = LEFT_MARGIN + self.start_pos
        y = self.gene_number * SPACE_BETWEEN_GENES + SPACE_BTW_HEADER_AND_GENE + TOP_MARGIN
        context.set_line_width(20)
        context.set_source_rgba(self.color[0][0], self.color[0][1], self.color[0][2], self.color[0][3])
        context.move_to(x ,y)
        context.line_to(x+self.length, y)
        context.stroke()

class GeneGroup:
    '''Parent class that controls the drawing of all pieces for a gene'''
    def __init__(self, gene_number):
        self.gene_number = gene_number
        self.motifs = []
        
    def draw(self): 
        self.header.draw()
        self.gene.draw()
        self.exon.draw()
        for motifobj_list in self.motifs:
            for motif in motifobj_list:
                motif.draw()



def generate_gene_groups_inator():
    '''generates GeneGroup objects with gene_numbers starting at 1'''
    gene_groups = list()
    for i in range(len(GENE_SEQUENCE_DICT.keys())): #the keys are the headers
        gene_groups.append(GeneGroup(i))
    return(gene_groups)

GENE_GROUPS = generate_gene_groups_inator()

def populate_objects_inator():
    '''populate GeneGroups and all subclassess'''
    for gene_group in GENE_GROUPS:
        gene_group.header = Header(list(GENE_SEQUENCE_DICT.keys())[gene_group.gene_number], gene_group.gene_number) #finding header

        gene_group.gene = Gene(len(list(GENE_SEQUENCE_DICT.values())[gene_group.gene_number]), list(GENE_SEQUENCE_DICT.values())[gene_group.gene_number], gene_group.gene_number) #finding gene

        exon_start_location = re.search('[A-Z]+', gene_group.gene.sequence).start() +1
        exon_end_location = re.search('[A-Z]+', gene_group.gene.sequence).end()
        exon_length = exon_end_location - exon_start_location
        gene_group.exon = Exon(exon_start_location, exon_length, gene_group.gene_number)

        motif_counter = 0
        for motif in MOTIF_REGEX_DICT.values():
            motif_objects = []
            motif_counter += 1
            motif_color = [(COLORS_INTS_LIST[motif_counter-1][0]/255, COLORS_INTS_LIST[motif_counter-1][1]/255, COLORS_INTS_LIST[motif_counter-1][2]/255, 0.7)]
            motif_starts_list, motif_lengths_list = motif_coordinates_inator(motif, gene_group.gene.sequence)
            for i in range(len(motif_starts_list)):
                motif_objects.append(Motif(motif_starts_list[i], motif_lengths_list[i], motif_color, gene_group.gene_number))
            gene_group.motifs.append(motif_objects)
        #Marking the motifs

populate_objects_inator()

### Getting shape of data (seq lengths, and # of genes)
seq_lengths = list()
for seq in (GENE_SEQUENCE_DICT.values()):
    seq_lengths.append(len(seq))
max_seq_length = max(seq_lengths)
number_of_genes = len(GENE_SEQUENCE_DICT.keys())

### initialize pycairo 'context' object (the surface we draw on)
figure_surface = cairo.SVGSurface(OUTPUT_FILE_NAME, max_seq_length + 100, number_of_genes * 75 + 100)
context = cairo.Context(figure_surface)
context.set_source_rgba(0, 0, 0, 1)
context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL)
context.set_font_size(12)

def pretty_picture_painter():
    '''populate GeneGroups and all subclassess'''
    for gene_group in GENE_GROUPS:
        gene_group.draw()
    #draws most of the figure

    context.set_source_rgba(0,0,0,1)
    x = LEFT_MARGIN
    y =len(GENE_SEQUENCE_DICT.values()) * (SPACE_BETWEEN_GENES +1) +20
    context.move_to(x,y)
    context.set_font_size(8)
    context.show_text('Legend')
    #Writing the word legend

    legend_counter = 0
    for motif in MOTIF_REGEX_DICT.keys():
        legend_counter += 1
        x = LEFT_MARGIN
        y = len(GENE_SEQUENCE_DICT.values()) * (SPACE_BETWEEN_GENES +1) + legend_counter * SPACE_BTW_LEGEND_MOTIFS + 20
        context.move_to(x,y)
        context.set_line_width(12)
        context.set_source_rgba(COLORS_INTS_LIST[legend_counter-1][0]/255,COLORS_INTS_LIST[legend_counter-1][1]/255,COLORS_INTS_LIST[legend_counter-1][2]/255,1)
        context.line_to(x+5, y)
        context.stroke()
        context.move_to(x+10, y)
        context.set_source_rgba(0,0,0,1)
        context.show_text(motif)
        context.stroke()

pretty_picture_painter()