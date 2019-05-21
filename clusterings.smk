_data_folder_ = "data/exp{exp}/ground"
_out_folder_ = "data/exp{exp}"

_sim_ = 50

_clusters_ = 100

_exp_ = [1,2,3]

tools = {
    'affinity': 'clustering_methods/affinity.py',
    'agglomerative': 'clustering_methods/agglomerative.py',
    'birch': 'clustering_methods/birch.py',
    'kmeans': 'clustering_methods/kmeans.py',
    'spectral': 'clustering_methods/spectral.py',
    'celluloid': 'clustering_methods/celluloid_test.py',
    'kmodes': 'clustering_methods/kmodes.py'
}

rule draw_plots:
    input:
        expand('plots/exp{exp}_plot_pr.pdf', exp=_exp_),

rule draw_accuracy:
    output:
        plot = 'plots/exp{exp}_plot_pr.pdf'
    input:
        flag = 'clustering.done',
        prgm = 'clustering_accuracy/clust_accuracies.py',
    shell:
        '''
            python3 {input.prgm} {wildcards.exp}
        '''

rule run_tools:
    output:
        touch('clustering.done')
    input:
        expand(_out_folder_ + '/affinity/sim_{sim}_scs_affinity_clusters.txt', sim = range(1, _sim_ + 1),
                exp=_exp_),
        expand(_out_folder_ + '/agglomerative/sim_{sim}_scs_agglomerative_clusters.txt', sim = range(1, _sim_ + 1),
                exp=_exp_),
        expand(_out_folder_ + '/birch/sim_{sim}_scs_birch_clusters.txt', sim = range(1, _sim_ + 1),
                exp=_exp_),
        expand(_out_folder_ + '/kmeans/sim_{sim}_scs_kmeans_clusters.txt', sim = range(1, _sim_ + 1),
                exp=_exp_),
        expand(_out_folder_ + '/spectral/sim_{sim}_scs_spectral_clusters.txt', sim = range(1, _sim_ + 1),
                exp=_exp_),
        expand(_out_folder_ + '/celluloid/sim_{sim}_scs_celluloid_clusters.txt', sim = range(1, _sim_ + 1),
                exp=_exp_),
        expand(_out_folder_ + '/kmodes/sim_{sim}_scs_kmodes_clusters.txt', sim = range(1, _sim_ + 1),
                exp=_exp_),


# ---------------------------------
# -------- Run Clusterings --------
# ---------------------------------

rule run_affinity:
    output:
        out_cluster = _out_folder_ + '/affinity/sim_{sim}_scs_affinity_clusters.txt',
        out_matrix = _out_folder_ + '/affinity/sim_{sim}_scs_affinity.matrix',
        out_mutations = _out_folder_ + '/affinity/sim_{sim}_scs_affinity.mutations',
    input:
        prgm = tools['affinity'],
        sim_file = _data_folder_ + '/sim_{sim}_scs.txt'
    params:
        k = _clusters_,
        out_folder = _out_folder_ + '/affinity'
    shell:
        '''
            python3 {input.prgm} -f {input.sim_file} \
                -k {params.k} \
                -o {params.out_folder} \
                -c mutations
        '''

rule run_agglomerative:
    output:
        out_cluster = _out_folder_ + '/agglomerative/sim_{sim}_scs_agglomerative_clusters.txt',
        out_matrix = _out_folder_ + '/agglomerative/sim_{sim}_scs_agglomerative.matrix',
        out_mutations = _out_folder_ + '/agglomerative/sim_{sim}_scs_agglomerative.mutations',
    input:
        prgm = tools['agglomerative'],
        sim_file = _data_folder_ + '/sim_{sim}_scs.txt'
    params:
        k = _clusters_,
        out_folder = _out_folder_ + '/agglomerative'
    shell:
        '''
            python3 {input.prgm} -f {input.sim_file} \
                -k {params.k} \
                -o {params.out_folder} \
                -c mutations
        '''

rule run_birch:
    output:
        out_cluster = _out_folder_ + '/birch/sim_{sim}_scs_birch_clusters.txt',
        out_matrix = _out_folder_ + '/birch/sim_{sim}_scs_birch.matrix',
        out_mutations = _out_folder_ + '/birch/sim_{sim}_scs_birch.mutations',
    input:
        prgm = tools['birch'],
        sim_file = _data_folder_ + '/sim_{sim}_scs.txt'
    params:
        k = _clusters_,
        out_folder = _out_folder_ + '/birch'
    shell:
        '''
            python3 {input.prgm} -f {input.sim_file} \
                -k {params.k} \
                -o {params.out_folder} \
                -c mutations
        '''

rule run_kmeans:
    output:
        out_cluster = _out_folder_ + '/kmeans/sim_{sim}_scs_kmeans_clusters.txt',
        out_matrix = _out_folder_ + '/kmeans/sim_{sim}_scs_kmeans.matrix',
        out_mutations = _out_folder_ + '/kmeans/sim_{sim}_scs_kmeans.mutations',
    input:
        prgm = tools['kmeans'],
        sim_file = _data_folder_ + '/sim_{sim}_scs.txt'
    params:
        k = _clusters_,
        out_folder = _out_folder_ + '/kmeans'
    shell:
        '''
            python3 {input.prgm} -f {input.sim_file} \
                -k {params.k} \
                -o {params.out_folder} \
                -c mutations
        '''

rule run_spectral:
    output:
        out_cluster = _out_folder_ + '/spectral/sim_{sim}_scs_spectral_clusters.txt',
        out_matrix = _out_folder_ + '/spectral/sim_{sim}_scs_spectral.matrix',
        out_mutations = _out_folder_ + '/spectral/sim_{sim}_scs_spectral.mutations',
    input:
        prgm = tools['spectral'],
        sim_file = _data_folder_ + '/sim_{sim}_scs.txt'
    params:
        k = _clusters_,
        out_folder = _out_folder_ + '/spectral'
    shell:
        '''
            python3 {input.prgm} -f {input.sim_file} \
                -k {params.k} \
                -o {params.out_folder} \
                -c mutations
        '''

rule run_celluloid:
    output:
        out_cluster = _out_folder_ + '/celluloid/sim_{sim}_scs_celluloid_clusters.txt',
        out_matrix = _out_folder_ + '/celluloid/sim_{sim}_scs_celluloid.matrix',
        out_mutations = _out_folder_ + '/celluloid/sim_{sim}_scs_celluloid.mutations',
    input:
        flag = 'kmodes/README.rst',
        prgm = tools['celluloid'],
        sim_file = _data_folder_ + '/sim_{sim}_scs.txt'
    params:
        k = _clusters_,
        out_folder = _out_folder_ + '/celluloid'
    shell:
        '''
            python3 {input.prgm} -f {input.sim_file} \
                -k {params.k} \
                -o {params.out_folder} \
                -c mutations
        '''

rule run_kmodes:
    output:
        out_cluster = _out_folder_ + '/kmodes/sim_{sim}_scs_kmodes_clusters.txt',
        out_matrix = _out_folder_ + '/kmodes/sim_{sim}_scs_kmodes.matrix',
        out_mutations = _out_folder_ + '/kmodes/sim_{sim}_scs_kmodes.mutations',
    input:
        flag = 'kmodes/README.rst',
        prgm = tools['kmodes'],
        sim_file = _data_folder_ + '/sim_{sim}_scs.txt'
    params:
        k = _clusters_,
        out_folder = _out_folder_ + '/kmodes'
    shell:
        '''
            python3 {input.prgm} -f {input.sim_file} \
                -k {params.k} \
                -o {params.out_folder} \
                -c mutations
        '''

# obtain kmodes from git repo
rule obtain_kmodes :
    output : 'kmodes/README.rst'
    shell : 'git clone https://github.com/nicodv/kmodes.git'
