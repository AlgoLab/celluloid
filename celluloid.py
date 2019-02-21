import sys
import argparse
import numpy as np
from kmodes.kmodes.kmodes import KModes
from kmodes.kmodes.util.dissim import matching_dissim, ng_dissim
from collections import defaultdict

# a 'weak' form of matching for (e.g., SCS) matrices --- in {0,1,2}
# --- which is the same as the matching_dissim, but does not count
# missing values --- 2's in this case
def conflict_dissim(a, b, **_) :

    v = np.vectorize(lambda ai, bi : ai != 2 and bi != 2 and ai != bi)

    return np.sum(v(a,b), axis = 1)

# likelihood that column a orignated from centroid b, according to
# false negative/positive rates, same as I from E in SASC, expressed
# as a dissimilarity
alpha = 0.1
beta = 0.001 # default values
def likelihood_dissim(a, b, **_) :

    # c.f. P(I|E) from SASC
    def p(i,e) :

        if i == 2 : return 1. # missing value, no contribution

        if e == 0 :
            return 1. - beta if i == 0 else beta
        if e == 1 :
            return alpha if i == 0 else 1. - alpha

        # case not in SASC, but we need to handle it here, since e can
        # be 2, we want a very lowwww.. probability when i is not 2
        return min(alpha, beta) / 1000. # e == 2 and i != 2

    # note: since we want a _dissimilarity_, we want this to be low
    # when the (log) likelihood is high, and vice versa
    v = np.vectorize(lambda i,e : -1. * np.log(p(i,e)))

    return np.sum(v(a,b), axis = 1)

# choose dissimilarity function
def choose_dissim(dstr) :
    if dstr == 'matching' :
        return matching_dissim
    if dstr == 'ng' :
        return ng_dissim
    if dstr == 'conflict' :
        return conflict_dissim
    if dstr == 'likelihood' :
        return likelihood_dissim

    return None

# gather false postive/negative rates from a file
def obtain_rates(filename) :

    rates = []
    for line in open(filename, 'r') :
        if 'rate:' in line :
            rates.append(float(line.split()[-1]))

    a,b,c = rates

    return a,b

#
# Parser
#----------------------------------------------------------------------

parser = argparse.ArgumentParser(description = '''

   run k-modes (https://github.com/nicodv/kmodes) on the columns of a
   datset

''')

parser.add_argument(
    '-i', '--input',
    metavar = 'INPUT', dest = 'input',
    type = str, default = sys.stdin,
    help = 'input file')
parser.add_argument(
    '-m', '--method',
    metavar = 'METHOD', dest = 'method',
    type = str, default = 'cao',
    help = 'initialization method')
parser.add_argument(
    '-d', '--dissim',
    metavar = 'DISSIM', dest = 'dissim',
    type = str, default = 'matching',
    help = 'dissimilarity function')
parser.add_argument(
    '-n', '--n_inits',
    metavar = 'N', dest = 'n',
    type = int, default = 10,
    help = 'number of iterations')
parser.add_argument(
    '-k', '--kmodes',
    metavar = 'K', dest = 'k',
    type = int, default = 10,
    help = 'number of modes')
parser.add_argument(
    '-r', '--rates',
    metavar = 'RATES', dest = 'rates',
    type = str,
    help = 'file containing false negative/positive rates')
parser.add_argument(
    '-c', '--centroids',
    metavar = 'CENTROIDS', dest = 'centroids',
    type = str,
    help = 'file to which the centroids are output')
parser.add_argument(
    '-v', '--verbose',
    dest = 'verbose',
    action = 'store_true',
    help = 'verbose')

args = parser.parse_args()

#
# Main
#----------------------------------------------------------------------

# false negative/positve rates
if args.dissim == 'likelihood' :
    alpha, beta = obtain_rates(args.rates)
    print('FN rate (alpha) =', alpha)
    print('FP rate (beta) =', beta)

# load (columns of) dataset
a = np.loadtxt(args.input, dtype = 'int').transpose()

# to trigger the setting of n_init to 1 (small bug, really)
method = args.method
if args.method.lower() == 'cao' :
	method = 'Cao'

km = KModes(
    n_clusters = args.k,
    cat_dissim = choose_dissim(args.dissim),
    init = method,
    n_init = args.n,
    verbose = 1 if args.verbose else 0)

clusters = km.fit_predict(a)

# gather clusters
labels = km.labels_
d = defaultdict(list)
i = 1
for label in labels :
    d[int(label)].append(str(i))
    i += 1

# output clusters in the agreed format
for key in d :
    print(key, '"'+','.join(d[key])+'"', sep = '\t', file = sys.stderr)

# Print the cluster centroids
c_file = open(args.centroids, 'w')
cs = km.cluster_centroids_
for c in cs :
    print(*c, file = c_file)
c_file.close()
