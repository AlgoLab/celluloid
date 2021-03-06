import numpy as np
from sklearn.cluster import AgglomerativeClustering

import argparse, os, sys, errno

parser = argparse.ArgumentParser(description='Agglomerative clustering')
parser.add_argument('-f', '--file', type=str, required=True,
                    help='SCS matrix')
parser.add_argument('-k', type=int, required=True,
                    help='number of features')
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

def most_common(lst):
    return max(set(lst), key=lst.count)

def cluster_and_output(k, matrix, clust_type, inputpath, outdir):
    model = AgglomerativeClustering(n_clusters=k, affinity="manhattan", 
        linkage="average").fit(matrix)
    # print(matrix.shape)

    from collections import defaultdict
    labels_to_index = defaultdict(list)

    # print(model.labels_)
    print(len(model.labels_))

    for ix, label in enumerate(model.labels_):
        labels_to_index[label].append(ix)
    
    print(len(labels_to_index))

    filename = os.path.splitext(os.path.basename(inputpath))[0]
    outfile = os.path.join(outdir, filename)

    with open('{0}_agglomerative.mutations'.format(outfile), 'w+') as file_out:
            for cluster in sorted(labels_to_index):
                file_out.write('{}\n'.format(
                    ','.join([ str(x+1) for x in labels_to_index[cluster]])
                ))
    
    with open('{0}_agglomerative_clusters.txt'.format(outfile), 'w+') as file_out:
        for cluster in sorted(labels_to_index):
            file_out.write('{0}\t"{1}"\n'.format(
                cluster, ','.join([ str(x+1) for x in labels_to_index[cluster]])
            ))

    out_matrix = []
    for cluster in sorted(labels_to_index):
        vectors = list(matrix[ix,:] for ix in labels_to_index[cluster])
        str_v = [''.join(list(map(str, list(map(int,v.tolist()))))) for v in vectors]
        out_matrix.append(list(map(int, list(most_common(str_v)))))
    
    out_matrix = np.transpose(np.array(out_matrix))
    np.savetxt('{0}_agglomerative.matrix'.format(outfile), out_matrix, fmt='%d', delimiter=' ')


    # print(len(labels_to_index))
    # from collections import defaultdict
    # stats = defaultdict(int)

    # print(out_matrix.shape)
    # for i in range(out_matrix.shape[0]):
    #     for j in range(out_matrix.shape[1]):
    #         stats[out_matrix[i,j]] += 1
    # print(stats)

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