
# SETTINGS TO BE PUT HERE

_exp_ = [1,2,3]

_clustering_folder_ = 'data/exp{exp}/{cluster}'

_output_folder_ = 'data/inference/exp{exp}'

_clusters_ = ['affinity', 'birch', 'agglomerative', 'spectral', 'kmeans', 'nocl',
                'celluloid', 'kmodes']
_sim_ = 50

convert = {
    'to_scite': 'inference_tools/convert/to_scite.py',
    'to_sphyr': 'inference_tools/convert/to_sphyr.py'
}

tools = {
    'scite': 'inference_tools/tools/scite',
    'sphyr_run': 'inference_tools/tools/sphyr/kDPFC',
    'sphyr_vis': 'inference_tools/tools/sphyr/visualize'
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
