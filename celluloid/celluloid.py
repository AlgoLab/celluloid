from __future__ import generator_stop

import sys
import argparse
import numpy as np
from kmodes.kmodes import KModes
from collections import defaultdict

#
# Converter
#-----------------------------------------------------
from tatsu.exceptions import FailedParse
from tatsu.buffering import Buffer
from tatsu.parsing import Parser
from tatsu.parsing import tatsumasu, leftrec, nomemo
from tatsu.parsing import leftrec, nomemo  # noqa
from tatsu.util import re, generic_main  # noqa


KEYWORDS = {}  # type: ignore


class SASCBuffer(Buffer):
    def __init__(
        self,
        text,
        whitespace=re.compile('(?s)[ \\t\\r\\f\\v]+'),
        nameguard=None,
        comments_re=None,
        eol_comments_re=None,
        ignorecase=None,
        namechars='',
        **kwargs
    ):
        super().__init__(
            text,
            whitespace=whitespace,
            nameguard=nameguard,
            comments_re=comments_re,
            eol_comments_re=eol_comments_re,
            ignorecase=ignorecase,
            namechars=namechars,
            **kwargs
        )


class SASCParser(Parser):
    def __init__(
        self,
        whitespace=re.compile('(?s)[ \\t\\r\\f\\v]+'),
        nameguard=None,
        comments_re=None,
        eol_comments_re=None,
        ignorecase=None,
        left_recursion=True,
        parseinfo=True,
        keywords=None,
        namechars='',
        tokenizercls=SASCBuffer,
        **kwargs
    ):
        if keywords is None:
            keywords = KEYWORDS
        super().__init__(
            whitespace=whitespace,
            nameguard=nameguard,
            comments_re=comments_re,
            eol_comments_re=eol_comments_re,
            ignorecase=ignorecase,
            left_recursion=left_recursion,
            parseinfo=parseinfo,
            keywords=keywords,
            namechars=namechars,
            tokenizercls=tokenizercls,
            **kwargs
        )

    @tatsumasu()
    def _start_(self):  # noqa
        self._file_()
        self._check_eof()

    @tatsumasu()
    def _cell_(self):  # noqa
        with self._choice():
            with self._option():
                self._token('0')
                self._cut()
            with self._option():
                self._token('1')
                self._cut()
            with self._option():
                self._token('2')
                self._cut()
            self._error('expecting one of: 0 1 2')

    @tatsumasu()
    def _row_(self):  # noqa

        def block0():
            self._cell_()
        self._closure(block0)

    @tatsumasu()
    def _file_(self):  # noqa

        def sep0():
            with self._group():
                self._token('\n')

        def block0():
            self._row_()
        self._gather(block0, sep0)


class SASCSemantics(object):
    def start(self, ast):  # noqa
        return ast

    def cell(self, ast):  # noqa
        return ast

    def row(self, ast):  # noqa
        return ast

    def file(self, ast):  # noqa
        return ast

class SCITEBuffer(Buffer):
    def __init__(
        self,
        text,
        whitespace=re.compile('(?s)[ \\t\\r\\f\\v]+'),
        nameguard=None,
        comments_re=None,
        eol_comments_re=None,
        ignorecase=None,
        namechars='',
        **kwargs
    ):
        super().__init__(
            text,
            whitespace=whitespace,
            nameguard=nameguard,
            comments_re=comments_re,
            eol_comments_re=eol_comments_re,
            ignorecase=ignorecase,
            namechars=namechars,
            **kwargs
        )


class SCITEParser(Parser):
    def __init__(
        self,
        whitespace=re.compile('(?s)[ \\t\\r\\f\\v]+'),
        nameguard=None,
        comments_re=None,
        eol_comments_re=None,
        ignorecase=None,
        left_recursion=True,
        parseinfo=True,
        keywords=None,
        namechars='',
        tokenizercls=SCITEBuffer,
        **kwargs
    ):
        if keywords is None:
            keywords = KEYWORDS
        super().__init__(
            whitespace=whitespace,
            nameguard=nameguard,
            comments_re=comments_re,
            eol_comments_re=eol_comments_re,
            ignorecase=ignorecase,
            left_recursion=left_recursion,
            parseinfo=parseinfo,
            keywords=keywords,
            namechars=namechars,
            tokenizercls=tokenizercls,
            **kwargs
        )

    @tatsumasu()
    def _start_(self):  # noqa
        self._file_()
        self._check_eof()

    @tatsumasu()
    def _cell_(self):  # noqa
        with self._choice():
            with self._option():
                self._token('0')
                self._cut()
            with self._option():
                self._token('1')
                self._cut()
            with self._option():
                self._token('2')
                self._cut()
            with self._option():
                self._token('3')
                self._cut()
            self._error('expecting one of: 0 1 2 3')

    @tatsumasu()
    def _row_(self):  # noqa

        def block0():
            self._cell_()
        self._closure(block0)

    @tatsumasu()
    def _file_(self):  # noqa

        def sep0():
            with self._group():
                self._token('\n')

        def block0():
            self._row_()
        self._gather(block0, sep0)


class SCITESemantics(object):
    def start(self, ast):  # noqa
        return ast

    def cell(self, ast):  # noqa
        return ast

    def row(self, ast):  # noqa
        return ast

    def file(self, ast):  # noqa
        return ast

class SPHYRBuffer(Buffer):
    def __init__(
        self,
        text,
        whitespace=re.compile('(?s)[ \\t\\r\\f\\v]+'),
        nameguard=None,
        comments_re=None,
        eol_comments_re='#([^\\n]*?)$',
        ignorecase=None,
        namechars='',
        **kwargs
    ):
        super().__init__(
            text,
            whitespace=whitespace,
            nameguard=nameguard,
            comments_re=comments_re,
            eol_comments_re=eol_comments_re,
            ignorecase=ignorecase,
            namechars=namechars,
            **kwargs
        )


class SPHYRParser(Parser):
    def __init__(
        self,
        whitespace=re.compile('(?s)[ \\t\\r\\f\\v]+'),
        nameguard=None,
        comments_re=None,
        eol_comments_re='#([^\\n]*?)$',
        ignorecase=None,
        left_recursion=True,
        parseinfo=True,
        keywords=None,
        namechars='',
        tokenizercls=SPHYRBuffer,
        **kwargs
    ):
        if keywords is None:
            keywords = KEYWORDS
        super().__init__(
            whitespace=whitespace,
            nameguard=nameguard,
            comments_re=comments_re,
            eol_comments_re=eol_comments_re,
            ignorecase=ignorecase,
            left_recursion=left_recursion,
            parseinfo=parseinfo,
            keywords=keywords,
            namechars=namechars,
            tokenizercls=tokenizercls,
            **kwargs
        )

    @tatsumasu()
    def _start_(self):  # noqa
        self._cell_number_line_()
        self._token('\n')
        self._snv_number_line_()
        self._token('\n')
        self._file_()
        self._check_eof()

    @tatsumasu()
    def _cell_number_line_(self):  # noqa
        self._pattern('\\d+')

    @tatsumasu()
    def _snv_number_line_(self):  # noqa
        self._pattern('\\d+')

    @tatsumasu()
    def _cell_(self):  # noqa
        with self._choice():
            with self._option():
                self._token('0')
                self._cut()
            with self._option():
                self._token('1')
                self._cut()
            with self._option():
                self._token('-1')
                self._cut()
            self._error('expecting one of: -1 0 1')

    @tatsumasu()
    def _row_(self):  # noqa

        def block0():
            self._cell_()
        self._closure(block0)

    @tatsumasu()
    def _file_(self):  # noqa

        def sep0():
            with self._group():
                self._token('\n')

        def block0():
            self._row_()
        self._gather(block0, sep0)


class SPHYRSemantics(object):
    def start(self, ast):  # noqa
        return ast

    def cell_number_line(self, ast):  # noqa
        return ast

    def snv_number_line(self, ast):  # noqa
        return ast

    def cell(self, ast):  # noqa
        return ast

    def row(self, ast):  # noqa
        return ast

    def file(self, ast):  # noqa
        return ast

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


# the conflict dissimilarity measure
def conflict_dissim(a, b, **_) :
    v = np.vectorize(lambda ai, bi : ai != 2 and bi != 2 and ai != bi)
    return np.sum(v(a,b), axis = 1)

def main():

    #
    # Parser
    #----------------------------------------------------------------------

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(
        dest="command",
        title='subcommands',
        description='valid subcommands')

    #
    # Translator_subparser
    #----------------------------------------------------------------------

    convert_parser = subparsers.add_parser("convert", description = '''

        convert input files to different formats

    ''')

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

    convert_parser.add_argument(
        "file",
        help="Input file")

    #
    # Cluster_subparser
    #----------------------------------------------------------------------

    cluster_parser = subparsers.add_parser("cluster", description = '''

    celluloid: clustering single cell sequencing data around centroids

    ''')

    cluster_parser.add_argument(
        '-i', '--input',
        metavar = 'INPUT', dest = 'input',
        type = str, default = sys.stdin,
        help = 'input file',
        required=True)

    cluster_parser.add_argument(
        '-n', '--n_inits',
        metavar = 'N', dest = 'n',
        type = int, default = 10,
        help = 'number of runs')

    cluster_parser.add_argument(
        '-k', '--kmodes',
        metavar = 'K', dest = 'k',
        type = int, required = True,
        help = 'number of modes')

    cluster_parser.add_argument(
        '-l', '--labels',
        metavar = 'LABELS', dest = 'labels',
        type = str, default = None,
        help = 'label file')

    cluster_parser.add_argument('-o',
        '--outdir', action='store', 
        type=str, required=True,
        help='output directory.')

    cluster_parser.add_argument(
        '-v', '--verbose',
        dest = 'verbose',
        action = 'store_true',
        help = 'verbose')

    args = parser.parse_args()
    #
    # Main
    #----------------------------------------------------------------------

    if args.command == "convert":
        convert(args)
    elif args.command == 'cluster':
        import os, errno

        try:
            os.makedirs(args.outdir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(args.outdir):
                pass
            else:
                raise

        # load (columns of) dataset
        a = np.loadtxt(args.input, dtype='int').transpose()

        # load (or create) labels
        snv_labels = list()
        if args.labels:
            with open(args.labels, 'r') as fin:
                for line in fin:
                    snv_labels.append(line.strip())
            if len(snv_labels) != a.shape[0]:
                sys.exit('SCS file and LABELS do not have the same number of mutations.')
        else:
            snv_labels = [str(x) for x in range(1, a.shape[0] + 1)]

        # run k-modes
        km = KModes(
            n_clusters=args.k,
            cat_dissim=conflict_dissim,
            init='huang',
            n_init=args.n,
            verbose=1 if args.verbose else 0)

        clusters = km.fit_predict(a)

        # obtain the clusters (of labels)
        labels = km.labels_
        d = defaultdict(list)
        i = 0
        for label in labels:
            d[int(label)].append(snv_labels[i])
            i += 1

        # store (only the non-empty) cluster centroids
        cs = km.cluster_centroids_
        centroids = list()
        for key in d:
            centroids.append(cs[key])

        out_matrix = np.array(centroids).transpose()

        # output clustered matrix (the non-empty centroids)
        filename, ext = os.path.splitext(os.path.basename(args.input))
        scs_outfile = os.path.join(args.outdir, filename + '_clustered' + ext)
        np.savetxt(scs_outfile, out_matrix, fmt='%d', delimiter=' ')

        # output mutation names
        filename, ext = 'LABELS', '.txt'
        if args.labels:
            filename, ext = os.path.splitext(os.path.basename(args.labels))
        labels_outfile = os.path.join(args.outdir, filename + '_clustered' + ext)

        with open(labels_outfile, 'w+') as fout:
            for key in sorted(d):
                fout.write('%s\n' % ','.join(d[key]))
    else:
        parser.print_help()

if __name__ == '__main__':
    main()