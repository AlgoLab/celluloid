import numpy as np
from sklearn.cluster import KMeans

import argparse, os, sys, errno

parser = argparse.ArgumentParser(description='Single Cell Generator v2')
parser.add_argument('-f', '--file', type=str, required=True,
                    help='SCS matrix')
parser.add_argument('-k', type=int, required=True,
                    help='K value of k-means')
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
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(matrix)

    from collections import defaultdict
    cluster_groups = defaultdict(list)

    for j in range(matrix.shape[0]):
        cluster_groups[kmeans.labels_[j]].append(j)
    
    tot_rows = 0
    for cluster in cluster_groups:
        tot_rows += len(cluster_groups[cluster])

    filename = os.path.splitext(os.path.basename(inputpath))[0]
    outfile = os.path.join(outdir, filename)

    centroids  = kmeans.cluster_centers_
    out_matrix = list()
    for c in centroids:
        x = list(map(int, list(map(round, c))))
        out_matrix.append(x)

    out_matrix = np.transpose(np.array(out_matrix))
    
    np.savetxt('{0}_kmean.matrix'.format(outfile), out_matrix, fmt='%d', delimiter=' ')

    with open('{0}_kmean_clusters.txt'.format(outfile), 'w+') as file_out:
        for cluster in sorted(cluster_groups):
            file_out.write('{0}\t"{1}"\n'.format(
                cluster, ','.join([ str(x+1) for x in cluster_groups[cluster]])
            ))

    with open('{0}_kmean.mutations'.format(outfile), 'w+') as file_out:
        for cluster in sorted(cluster_groups):
            file_out.write('{0}\n'.format(
                ','.join([ str(x+1) for x in cluster_groups[cluster]])
            ))

    print('Done.')


if args.cluster == 'cells':
    print('Clustering cells')
    cluster_and_output(args.k, scs_matrix_input, args.cluster, args.file, args.outdir)
elif args.cluster == 'mutations':
    scs_matrix_input = np.transpose(scs_matrix_input)
    print('Clustering mutations')
    cluster_and_output(args.k, scs_matrix_input, args.cluster, args.file, args.outdir)
else:
    sys.exit('Something very wrong happened.')

