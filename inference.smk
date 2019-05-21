'''

IMPORTANT: this pipeline has been set up to run the downstream
inference methods, namely scite and sphyr, under the assumption that
they have been in installed in _toolsdir_ tools directory (below),
more specifically, that the associated excecutables have been
installed in tools dictionary of executables (below).  Now, provided
all dependencies for these tools have been installed correctly, this
pipeline will download, install and place these executables in the
right place, but this is only contingent on the proper setup of the
dependencies, which must be followed according to instructions at

    https://github.com/cbg-ethz/SCITE (for scite)

    https://github.com/elkebir-group/SPhyR (for sphyr)

only after all such dependencies are installed, should this pipeline
run rather automatically

'''

_toolsdir_ = 'inference_tools/tools' # tools directory

tools = { # executables
    'scite': _toolsdir_ + '/SCITE/scite',
    'sphyr_run': _toolsdir_ + '/SPhyR/build/kDPFC',
    'sphyr_vis': _toolsdir_ + '/SPhyR/build/visualize'
}

_exp_ = [1,2,3]

_clustering_folder_ = 'data/exp{exp}/{cluster}'

_output_folder_ = 'data/inference/exp{exp}'

_clusters_ = ['affinity birch agglomerative spectral kmeans nocl celluloid kmodes'.split()

_sim_ = 50

convert = {
    'to_scite': 'inference_tools/convert/to_scite.py',
    'to_sphyr': 'inference_tools/convert/to_sphyr.py'
}

_fn_ = 0.250000
_fp_ = 0.000010
_sphyr_k_ = 0

# ---------------------------------
# --------- Run tools -------------
# ---------------------------------
rule run_tools:
    input:
        expand(_output_folder_ + '/scite/{cluster}/sim_{sim}_scs_{cluster}_ml0.gv', 
            sim = range(1, _sim_ + 1), 
            cluster = _clusters_, 
            exp = _exp_),
        expand(_output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}.gv',
            sim = range(1, _sim_ + 1), 
            cluster = _clusters_, 
            exp = _exp_)


rule run_scite:
    output:
        scite_result = _output_folder_ + '/scite/{cluster}/sim_{sim}_scs_{cluster}_ml0.gv'
    input:
        prgm = tools['scite'],
        scite_scs = _output_folder_ + '/scite/{cluster}/sim_{sim}_scs_{cluster}.in',
        scite_names = _output_folder_ + '/scite/{cluster}/sim_{sim}_scs_{cluster}.geneNames'
    params:
        fn_rate = _fn_,
        fp_rate = _fp_
    benchmark:
        'benchmark/inference/exp{exp}/scite/{cluster}/sim_{sim}.benchmark.txt'
    run:
        with open(input.scite_scs) as mat:
            m = len(mat.readline().split())
            n = sum(1 for line in mat) + 1

        shell(
            '''
                {input.prgm} -i {input.scite_scs} -n {n} -m {m} \
                    -r 1 -l 900000 -fd {params.fp_rate} -ad {params.fn_rate} \
                        -names {input.scite_names} -max_treelist_size 1 -a
            ''')

rule run_sphyr:
    output:
        sphyr_vis = _output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}.gv',
        sphyr_out = _output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}.A',
    input:
        prgm_run = tools['sphyr_run'],
        prgm_vis = tools['sphyr_vis'],
        sphyr_scs = _output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}.input',
        sphyr_genenames = _output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}_SNV.labels',
        sphyr_cellnames = _output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}.cellNames',
    params:
        fn_rate = _fn_,
        fp_rate = _fp_,
        k = _sphyr_k_
    benchmark:
        'benchmark/inference/exp{exp}/sphyr/{cluster}/sim_{sim}.benchmark.txt' 
    shell:
        '''
            {input.prgm_run} {input.sphyr_scs} \
                -a {params.fp_rate} -b {params.fn_rate} \
                -k {params.k} > {output.sphyr_out} && \
            {input.prgm_vis} {output.sphyr_out} \
                -c {input.sphyr_genenames} \
                -t {input.sphyr_cellnames} > {output.sphyr_vis}
        '''

# ---------------------------------
# ------- Prepare sim data --------
# ---------------------------------

rule prepare_simulations:
    input:
        expand(_output_folder_ + '/scite/{cluster}/sim_{sim}_scs_{cluster}.in', 
            sim = range(1, _sim_ + 1), 
            cluster = _clusters_, 
            exp = _exp_),
        expand(_output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}.input', 
            sim = range(1, _sim_ + 1), 
            cluster = _clusters_, 
            exp = _exp_),
        expand(_output_folder_ + '/sasc/{cluster}/sim_{sim}_scs_{cluster}.matrix', 
            sim = range(1, _sim_ + 1), 
            cluster = _clusters_, 
            exp = _exp_),


# ---------------------------------
# ------- Convert Tools -----------
# ---------------------------------

rule convert_to_scite:
    output:
        scite_scs = _output_folder_ + '/scite/{cluster}/sim_{sim}_scs_{cluster}.in',
        scite_names = _output_folder_ + '/scite/{cluster}/sim_{sim}_scs_{cluster}.geneNames'
    input:
        prgm = convert['to_scite'],
        sasc_scs = _clustering_folder_ + '/sim_{sim}_scs_{cluster}.matrix',
        mutations = _clustering_folder_ + '/sim_{sim}_scs_{cluster}.mutations',
    params:
        outdir = _output_folder_ + '/scite/{cluster}'
    shell:
        '''
            python3 {input.prgm} -f {input.sasc_scs} -o {params.outdir} && \
            cp {input.mutations} {params.outdir}
        '''

rule convert_to_sphyr:
    output:
        sphyr_scs = _output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}.input',
        sphyr_genenames = _output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}_SNV.labels',
        sphyr_cellnames = _output_folder_ + '/sphyr/{cluster}/sim_{sim}_scs_{cluster}.cellNames',
    input:
        prgm = convert['to_sphyr'],
        sasc_scs = _clustering_folder_ + '/sim_{sim}_scs_{cluster}.matrix',
        mutations = _clustering_folder_ + '/sim_{sim}_scs_{cluster}.mutations',
    params:
        outdir = _output_folder_ + '/sphyr/{cluster}'
    shell:
        '''
            python3 {input.prgm} -f {input.sasc_scs} -o {params.outdir} && \
            cp {input.mutations} {params.outdir}
        '''


#
# obtain and build scite
#----------------------------------------------------------------------
rule build_scite :
    output : tools['scite']
    input : _toolsdir_ + '/SCITE/README.md'
    shell : 'cd {_toolsdir_}/SCITE && g++ -std=c++11 *.cpp -o scite'

# obtain scite
rule download_scite :
    output : _toolsdir_ + '/SCITE/README.md'
    message : 'downloading SCITE to {_toolsdir_}'
    shell : 'cd {_toolsdir_} && git clone https://github.com/cbg-ethz/SCITE'


# obtain and build sphyr
#----------------------------------------------------------------------
rule build_sphyr :
    output :
        run = tools['sphyr_run'],
        vis = tools['sphyr_vis']

    input : _toolsdir_ + '/SPhyR/README.md'
    shell : '''

   cd {_toolsdir_}/SPhyR && mkdir -p build && cd build && cmake .. && make '''

# obtain sphyr
rule download_sphyr :
    output : _toolsdir_ + '/SPhyR/README.md'
    message : 'downloading SPhyR to {_toolsdir_}'
    shell : 'cd {_toolsdir_} && git clone https://github.com/elkebir-group/SPhyR.git'
