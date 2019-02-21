import numpy as np
import sys

centroids = []
# read in the centroids file
for line in open(sys.argv[1], 'r') :
	centroids.append(line.split())

kept = []
# read in the clusters (modes) file
for line in open(sys.argv[2], 'r') :
	i = int(line.split()[0])

	# only keep the centroids with non-empty clusters
	kept.append(centroids[i])

a = np.array(kept, dtype = 'int').transpose()

np.savetxt(sys.argv[3], a, fmt = '%d', delimiter = ' ')
