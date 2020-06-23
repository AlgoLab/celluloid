import sys
import argparse
import numpy as np
from kmodes.kmodes import KModes
from collections import defaultdict
from inference_tools.convert import Translator

# the conflict dissimilarity measure
def conflict_dissim(a, b, **_) :
    v = np.vectorize(lambda ai, bi : ai != 2 and bi != 2 and ai != bi)
    return np.sum(v(a,b), axis = 1)

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
    help = 'number of iterations')

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
    Translator.convert(args)
else:
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
