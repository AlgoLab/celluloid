import sys
sys.path.insert(0, 'kmodes/')

import numpy as np
from kmodes.kmodes import KModes

# the conflict dissimilarity measure
def conflict_dissim(a, b, **_) :
    v = np.vectorize(lambda ai, bi : ai != 2 and bi != 2 and ai != bi)
    return np.sum(v(a,b), axis = 1)

import argparse, os, sys, errno

parser = argparse.ArgumentParser(description='K-means clustering')
parser.add_argument('-f', '--file', type=str, required=True,
                    help='SCS matrix')
parser.add_argument('-k', type=int, required=True,
                    help='K value of celluloid')
parser.add_argument('-n', type=int, default=10,
                    help='n_init')
parser.add_argument('-c', '--cluster', required=True,
                    choices=['cells', 'mutations', 'both'],
                    help='Cluster either cells or mutations')
parser.add_argument('-o', '--outdir', type=str, required=True,
                    help='output path')

args = parser.parse_args()

scs_matrix_input = np.loadtxt(args.file, dtype='d', delimiter=' ')

# create output directory if it doesn't already exists

try:
    os.makedirs(args.outdir)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(args.outdir):
        pass
    else:
        raise

def cluster_and_output(k, matrix, clust_type, inputpath, outdir):
    km = KModes(
    n_clusters = args.k,
    cat_dissim = conflict_dissim,
    init = 'huang',
    n_init = args.n,
    verbose = 1)
    km.fit(matrix)

    from collections import defaultdict
    cluster_groups = defaultdict(list)

    for j in range(matrix.shape[0]):
        cluster_groups[km.labels_[j]].append(j)
    
    tot_rows = 0
    for cluster in cluster_groups:
        tot_rows += len(cluster_groups[cluster])

    filename = os.path.splitext(os.path.basename(inputpath))[0]
    outfile = os.path.join(outdir, filename)

    centroids  = km.cluster_centroids_
    out_matrix = list()
    for ix_c, c in enumerate(centroids):
        if ix_c in cluster_groups:
            x = list(map(int, list(map(round, c))))
            out_matrix.append(x)

    out_matrix = np.transpose(np.array(out_matrix))

    print(out_matrix.shape)
    print(len(cluster_groups))
    np.savetxt('{}_celluloid.matrix'.format(outfile), out_matrix, fmt='%d', delimiter=' ')

    with open('{}_celluloid_clusters.txt'.format(outfile), 'w+') as file_out:
        for cluster in sorted(cluster_groups):
            file_out.write('{0}\t"{1}"\n'.format(
                cluster, ','.join([ str(x+1) for x in cluster_groups[cluster]])
            ))

    with open('{}_celluloid.mutations'.format(outfile), 'w+') as file_out:
        for cluster in sorted(cluster_groups):
            file_out.write('{0}\n'.format(
                ','.join([ str(x+1) for x in cluster_groups[cluster]])
            ))

    print('Done.')


if args.cluster == 'cells':
    # print('Clustering cells')
    # cluster_and_output(args.k, scs_matrix_input, args.cluster, args.file, args.outdir)
    print('not fully implemented yet')
elif args.cluster == 'mutations':
    scs_matrix_input = np.transpose(scs_matrix_input)
    print('Clustering mutations')
    cluster_and_output(args.k, scs_matrix_input, args.cluster, args.file, args.outdir)
elif args.cluster == 'both':
    print('not implemented yet')
else:
    sys.exit('Something very wrong happened.')
