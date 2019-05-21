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

inp = np.genfromtxt(args.file, delimiter = DELIMITER, dtype=int)

matrix_input = np.transpose(inp)

for i in range(matrix_input.shape[0]):
	for j in range(matrix_input.shape[1]):
		if matrix_input[i, j] == 2:
			matrix_input[i, j] = 3

filename = os.path.splitext(os.path.basename(args.file))[0]
outfile = os.path.join(args.out, filename)

np.savetxt('%s.in' % outfile, matrix_input, delimiter=' ', fmt='%i')

with open('%s.geneNames' % outfile, 'w+') as names:
    for name in range(1, matrix_input.shape[0] + 1):
        names.write('%s\n' % name)
