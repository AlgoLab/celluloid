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

n = inp.shape[0]
m = 1
if len(inp.shape) > 1 :
        m = inp.shape[1]
else :
        inp = np.reshape(inp, (n,1))

matrix_input = np.transpose(inp)

for i in range(m):
	for j in range(n):
		if matrix_input[i, j] == 2:
			matrix_input[i, j] = 3

filename = os.path.splitext(os.path.basename(args.file))[0]
outfile = os.path.join(args.out, filename)

np.savetxt('%s.in' % outfile, matrix_input, delimiter=' ', fmt='%i')

with open('%s.geneNames' % outfile, 'w+') as names:
    for name in range(1, m + 1):
        names.write('%s\n' % name)
