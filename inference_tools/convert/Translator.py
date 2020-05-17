import argparse
import sys
from inference_tools.convert.parserSASC import SASCParser
from inference_tools.convert.parserSCITE import SCITEParser
from inference_tools.convert.parserSPHYR import SPHYRParser

from tatsu.exceptions import FailedParse


# command line options
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outputFormat", help="Format to translate the input to (default is SASC)",
                    choices=["SASC", "SCITE", "SPHYR"])
parser.add_argument("-i", "--inputFormat", help="Format to translate the input from (default is SASC)",
                    choices=["SASC", "SCITE", "SPHYR"])
parser.add_argument("--outfile", help="Output file (default is stdout)")
parser.add_argument("file", type=str, help="Input file")
args = parser.parse_args()

# used to choose whether to transpose the matrix or not in the translation
transpose = {
    "SASC": False,
    "SCITE": True,
    "SPHYR": False
}
# the character used to say that there is no information in each format, None might not be needed
no_info = {
    "SASC": "2",
    "SCITE": "3",
    "SPHYR": "-1",
    None: "2"
}
# gets the input and output format from command line, defaults to SASC if none are given
format_to = args.outputFormat if args.outputFormat else "SASC"
format_from = args.inputFormat if args.inputFormat else "SASC"

# gets output and input file names from command line
output_file_name = args.outfile
input_file_name = args.file

# reads the input file
with open(input_file_name, "r") as input_file:
    file_str = input_file.read()


parser_dict = {
    "SASC": SASCParser,
    "SCITE": SCITEParser,
    "SPHYR": SPHYRParser
}
parser = parser_dict[format_from]()
try:
    ast = parser.parse(file_str, rule_name='start')
except FailedParse:
    sys.exit("wrong input format")

# handles the header of the SPHYR format
if format_from == "SPHYR":
    # cells and SNVs maybe not needed
    cells_number = ast[0]
    SNVs_number = ast[2]
    ast = ast[4]

# remove empty lines
ast = [row for row in ast if row != []]

if format_from == "SCITE":
    # SCITE sometimes uses both 1 and 2 for heterozygous and homozygous mutations
    ast = [[1 if a == 2 else a for a in row] for row in ast]
# changes the values of the "no info" character from the one of the input format to the one of the output format
ast_changed = [[a if a != no_info[format_from] else no_info[format_to] for a in row] for row in ast]

# transposes the matrix if needed
if transpose[format_to] != transpose[format_from]:
    # there could be a better way to do this
    ast_changed = list(map(list, zip(*ast_changed)))

# writes the translated file
# this is really ugly, any better way to do it?
if output_file_name is None:
    if format_to == "SPHYR":
        sys.stdout.write(str(len(ast_changed)) + " #cells\n")
        sys.stdout.write(str(len(ast_changed[1])) + " #SNVs\n")
    for row in ast_changed:
        for cell in row:
            sys.stdout.write(cell)
            sys.stdout.write(" ")
        sys.stdout.write("\n")
else:
    with open(output_file_name, "w") as file:
        if format_to == "SPHYR":
            file.write(str(len(ast_changed)) + " #cells\n")
            file.write(str(len(ast_changed[1])) + " #SNVs\n")
        for row in ast_changed:
            for cell in row:
                file.write(cell)
                file.write(" ")
            file.write("\n")
