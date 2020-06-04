import argparse
import sys
from inference_tools.convert.parserSASC import SASCParser
from inference_tools.convert.parserSCITE import SCITEParser
from inference_tools.convert.parserSPHYR import SPHYRParser
from tatsu.exceptions import FailedParse


def parse_file(file_name, file_format):
    # reads the input file
    with open(file_name, "r") as input_file:
        file_str = input_file.read()
    parser_dict = {
        "SASC": SASCParser,
        "SCITE": SCITEParser,
        "SPHYR": SPHYRParser
    }
    file_parser = parser_dict[file_format]()
    try:
        ast = file_parser.parse(file_str, rule_name='start')
    except FailedParse:
        print("wrong input format")
        raise
    # remove empty lines
    ast = [row for row in ast if row != []]
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
    # handles the header of the SPHYR format
    if format1 == "SPHYR":
        # cells and SNVs maybe not needed
        cells_number = input_ast[0]
        SNVs_number = input_ast[2]
        input_ast = input_ast[4]

    # Change 1s and 2s in SCITE to only 1s, since they both represent a mutation
    if format1 == "SCITE":
        input_ast = [[1 if a == 2 else a for a in row] for row in input_ast]

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
            sys.stdout.write(str(len(ast_translated[1])) + " #SNVs\n")
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
                file.write(str(len(ast_translated[1])) + " #SNVs\n")
            for row in ast_translated:
                for cell in row:
                    file.write(cell)
                    file.write(" ")
                file.write("\n")


if __name__ == "__main__":
    # command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outputFormat", help="Format to translate the input to (default is SASC)",
                        choices=["SASC", "SCITE", "SPHYR"])
    parser.add_argument("-i", "--inputFormat", help="Format to translate the input from (default is SASC)",
                        choices=["SASC", "SCITE", "SPHYR"])
    parser.add_argument("--outfile", help="Output file (default is stdout)")
    parser.add_argument("file", type=str, help="Input file")
    args = parser.parse_args()

    # gets the input and output format from command line, defaults to SASC if none are given
    format_to = args.outputFormat if args.outputFormat else "SASC"
    format_from = args.inputFormat if args.inputFormat else "SASC"
    # gets output and input file names from command line
    output_file_name = args.outfile
    input_file_name = args.file

    # reads the file and parses it
    ast = parse_file(input_file_name, format_from)
    print(ast)
    # translates the parsed file into the desired format
    ast_changed = translate(ast, format_from, format_to)
    print(ast_changed)
    # writes the file
    write_file(ast_changed, output_file_name, format_to)
    write_file(ast_changed, None, format_to)
    test_ast = [['2', '1', '2', '0'], ['2', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '2', '0']]
    print(test_ast)
    test_ast2 = translate(test_ast, format_from, format_to)
    print(test_ast2)



