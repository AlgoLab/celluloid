# celluloid
Celluloid: clustering single cell sequencing data around centroids

How to install
-----------------
Install the dependecies of the tool, with `pip3 install kmodes` and then just clone this repository:

```bash
git clone https://github.com/AlgoLab/celluloid.git
```

Input formats
-----------------
**Single Cell file**

The input file (specified by the `-i` parameter) is expected to be a ternary matrix file where the rows represent the cells and the columns the mutations. Each cell must be separated by a space. Each cell of the matrix can be:

| Value of cell | Meaning |
| ------------- | ------------- |
| I[i,j] = 0    | Mutation *j* is not observed in cell *i*  |
| I[i,j] = 1    | Mutation *j* is observed in cell *i*  |
| I[i,j] = 2    | There is no information for mutation *j* in cell *i*, i.e. low coverage  |

**Mutations file**

This optional file specifies the name of the mutations (parameter `-l`). Each mutation's name must be on a different line (separated by `\n`), and the names are assigned to columns from left to right in the input file. If this file is not provided, mutations are progressively named from `1` to the total number of mutations.

Output formats
-----------------
Celluloid will output one _Single Cell file_ and one _Mutations file_ in the same format as the input but with a number of mutations reduced by the clustering. Labels of mutation are concatenated in the out label file.


How to run
-----------------
```
usage: celluloid.py [-h] -i INPUT [-m {huang,random}] [-d {conflict,matching}]
                    [-n N] -k K [-l LABELS] -o OUTDIR [-v]

Celluloid: clustering single cell sequencing data around centroids

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file
  -m {huang,random}, --method {huang,random}
                        initialization method
  -d {conflict,matching}, --dissim {conflict,matching}
                        dissimilarity measure
  -n N, --n_inits N     number of iterations
  -k K, --kmodes K      number of modes
  -l LABELS, --labels LABELS
                        label file
  -o OUTDIR, --outdir OUTDIR
                        output directory.
  -v, --verbose         verbose
```

Replicating the experiments
-----------------

The preprint of an article outlining a study on Celluloid can be found
at [https://doi.org/10.1101/586545](https://doi.org/10.1101/586545)

A workflow for reproducing the experiments of this study is described
as a "Snakefile", which needs to be run with
[snakemake](http://snakemake.bitbucket.org).

**Overview**

The workflow is designed to be run from this directory.  In principle,
you only need to run `snakemake` here and then the rest is done
automatically, but all dependencies need to be installed first.  To
install dependencies, see [How to install](#how-to-install), and then
install [sasc](https://github.com/sciccolella/sasc):

    git clone https://github.com/sciccolella/sasc
    cd sasc && gcc -o sasc *.c -std=gnu99 -g -lm

**Run the workflow**

Start the workflow with

    nice snakemake -p -j 16

Adjust the 16, which is the number of cores you want to use.

**Results**

The resulting files for each dataset, _e.g._,
`data/exp1/sim_1_scs.txt` will be reported in a directory named
`data/exp1/sim_1.{parameters}.celluloid` (resp.,
`data/exp1/sim_1.{parameters}.kmeans`) for the clustering -- and the
downstream run of sasc -- with celluloid (resp., kmeans), where
`{parameters}` are the parameters, _e.g._, initialization method,
dissimilarity measure, k, associated with corresponding clustering
method
