import argparse
import sys
from .parserSASC import SASCParser
from .parserSCITE import SCITEParser
from .parserSPHYR import SPHYRParser
from tatsu.exceptions import FailedParse


class NotAMatrix(Exception):
    pass


def parse_string(file_as_string, file_format):
    parser_dict = {
        "SASC": SASCParser,
        "SCITE": SCITEParser,
        "SPHYR": SPHYRParser
    }
    file_parser = parser_dict[file_format]()
    try:
        ast = file_parser.parse(file_as_string, rule_name='start')
    except FailedParse:
        raise
    if file_format == "SPHYR":
        ast = ast[4]
    # remove empty lines
    ast = [row for row in ast if row != []]
    for row in ast:
        if len(row) != len(ast[0]):
            raise NotAMatrix("Number of cells per row varies")
    return ast


def translate(input_ast, format1, format2):
    # the character used to say that there is no information in each format
    no_info = {
        "SASC": "2",
        "SCITE": "3",
        "SPHYR": "-1",
    }
    # used to choose whether to transpose the matrix or not in the translation
    transpose = {
        "SASC": False,
        "SCITE": True,
        "SPHYR": False
    }

    # Change 1s and 2s in SCITE to only 1s, since they both represent a mutation
    if format1 == "SCITE":
        input_ast = [['1' if a == '2' else a for a in row] for row in input_ast]

    # changes the values of the "no info" character from the one of the input format to the one of the output format
    ast_translated = [[a if a != no_info[format1] else no_info[format2] for a in row] for row in input_ast]

    # transposes the matrix if needed
    if transpose[format2] != transpose[format1]:
        ast_translated = list(map(list, zip(*ast_translated)))

    return ast_translated


def write_file(ast_translated, file_name, file_format):
    # writes the translated file
    if file_name is None:
        # SPHYR needs its header
        if file_format == "SPHYR":
            sys.stdout.write(str(len(ast_translated)) + " #cells\n")
            sys.stdout.write(str(len(ast_translated[0])) + " #SNVs\n")
        for row in ast_translated:
            for cell in row:
                sys.stdout.write(cell)
                sys.stdout.write(" ")
            sys.stdout.write("\n")
    else:
        with open(file_name, "w") as file:
            # SPHYR needs its header
            if file_format == "SPHYR":
                file.write(str(len(ast_translated)) + " #cells\n")
                file.write(str(len(ast_translated[0])) + " #SNVs\n")
            for row in ast_translated:
                for cell in row:
                    file.write(cell)
                    file.write(" ")
                file.write("\n")


def convert(arguments):
    # gets the input and output format from command line, defaults to SASC if none are given
    format_to = arguments.outputFormat if arguments.outputFormat else "SASC"
    format_from = arguments.inputFormat if arguments.inputFormat else "SASC"
    # gets output and input file names from command line
    output_file_name = arguments.outfile
    input_file_name = arguments.file

    # reads the input file
    with open(input_file_name, "r") as input_file:
        file_str = input_file.read()
    # parses the input file
    try:
        ast = parse_string(file_str, format_from)
    except FailedParse:
        print("wrong input format")
        sys.exit(0)
    except NotAMatrix:
        print("input is not a matrix")
        sys.exit(0)
    # translates the parsed file into the desired format
    ast_changed = translate(ast, format_from, format_to)
    # writes the file
    write_file(ast_changed, output_file_name, format_to)


if __name__ == "__main__":
    # command line options
    convert_parser = argparse.ArgumentParser()

    convert_parser.add_argument(
        "-o",
        "--outputFormat",
        help="Format to translate the input to (default is SASC)",
        choices=["SASC", "SCITE", "SPHYR"])

    convert_parser.add_argument(
        "-i",
        "--inputFormat",
        help="Format to translate the input from (default is SASC)",
        choices=["SASC", "SCITE", "SPHYR"])

    convert_parser.add_argument(
        "--outfile",
        help="Output file (default is stdout)")

    convert_parser.add_argument("file",
                                help="Input file")

    args = convert_parser.parse_args()
    convert(args)





