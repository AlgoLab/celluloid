#!/usr/bin/env python3
import numpy as np
import argparse, os

# ======== INFO
'''
Convert SASC input to SCITE
'''

# ======== COMMAND LINE THINGS
parser = argparse.ArgumentParser(description='converter', add_help=True)
parser.add_argument('-f', '--file', required = True,
					type = str,
					help = 'Input file')
parser.add_argument('-o', '--out', required = True,
					type = str,
					help = 'Out Dir')
args = parser.parse_args()

# ======== INPUT PROCESSING FOR REAL DATA [CSV]:
DELIMITER = ' '

matrix_input = np.genfromtxt(args.file, delimiter = DELIMITER, dtype=int)

MUTATIONS = matrix_input.shape[1]

for i in range(matrix_input.shape[0]):
	for j in range(matrix_input.shape[1]):
		if matrix_input[i, j] == 2:
			matrix_input[i, j] = -1

filename = os.path.splitext(os.path.basename(args.file))[0]
outfile = os.path.join(args.out, filename)

with open('%s.input' % outfile, 'w+') as fout:
    fout.write('%d #cells\n%d #SNVs\n' % (matrix_input.shape[0], MUTATIONS))
    for i in range(matrix_input.shape[0]):
        fout.write(' '.join(map(str, matrix_input[i,:])))
        fout.write('\n')


with open('%s_SNV.labels' % outfile, 'w+') as names:
    for name in range(1, MUTATIONS + 1):
        names.write('%s\n' % name)

with open('%s.cellNames' % outfile, 'w+') as names:
    for name in range(1, matrix_input.shape[0] + 1):
        names.write('cell%s\n' % name)