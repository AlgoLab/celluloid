import load_clusters
import sys
import numpy as np
import pandas as pd
import itertools
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import completeness_score
from sklearn.metrics.cluster import v_measure_score

def build_heatmap_matrix(labels, score):
    tot_clusters = len(labels)
    hm_matrix = np.zeros((tot_clusters, tot_clusters))
    
    first_columns = list()

    tot_simulations = len(labels[0])
    
    for i in range(hm_matrix.shape[0]):
        for j in range(i, hm_matrix.shape[1]):
        # for j in range(hm_matrix.shape[1]):
            values = list()
            for sim in range(tot_simulations):
                values.append(score(labels[i][sim][2], labels[j][sim][2]))
            hm_matrix[i][j] = np.average(values)
            if i == 0 and j != 0:
                first_columns.append(values)

    
    return hm_matrix, first_columns

def plot_heatmap(hm_matrix, axis_labels, title, outfile=None):
    import seaborn as sns
    import matplotlib.pyplot as plt

    hm_matrix = np.transpose(hm_matrix)

    mask = np.zeros_like(hm_matrix, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    for i in range(mask.shape[0]):
        mask[i,i] = False

    f, ax = plt.subplots(figsize=(7, 7))
    sns.heatmap(hm_matrix, cmap="Greens", mask=mask, vmin=0, center=0,
            square=True, 
            linewidths=.5, 
            # cbar_kws={"shrink": .5},
            annot=True, fmt='.2f',            
            xticklabels=axis_labels, yticklabels=axis_labels)
    
    plt.yticks(rotation=0)

    plt.title(title)

    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
    # plt.show()

def plot_heatmaps_compacted_2x2(hm_matrices, axis_labels, titles, outfile=None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # grid_kws = {"width_ratios": (.24, .24, .24, .24, .04)}
    f, ax = plt.subplots(2,2, figsize=(9, 9)) #, gridspec_kw=grid_kws)

    mask = np.zeros_like(hm_matrices[0], dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    for i in range(mask.shape[0]):
        mask[i,i] = False

    for ix, hm_matrix in enumerate(hm_matrices):

        hm_matrix = np.transpose(hm_matrix)
        
        if ix == 0 or ix == 2:
            yticks = axis_labels
        else:
            yticks = False
        # if ix == :
        #     cbar = False

        cbar = False
        cbar_ax = None

        if ix == 0: 
            a = 0
            b = 0
        if ix == 1:
            a = 0
            b = 1
        if ix == 2:
            a = 1
            b = 0
        if ix == 3:
            a = 1
            b = 1
        
        ax[a,b] = sns.heatmap(hm_matrix, cmap="Greens", mask=mask, vmin=0, center=0,
                ax=ax[a,b],
                square=True, 
                linewidths=.3, 
                cbar=False,
                annot=True, fmt='.2f', annot_kws={"size": 10},
                xticklabels=axis_labels, yticklabels=yticks)
        
        # ax[ix].set(title=titles[ix])
        ax[a,b].set_title(titles[ix], fontsize=12)

        plt.yticks(rotation=0)

    # plt.title(title)

    plt.tight_layout()
    # plt.show()
    if outfile:
        plt.savefig(outfile)

def plot_heatmaps_compacted(hm_matrices, axis_labels, titles, outfile=None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # grid_kws = {"width_ratios": (.24, .24, .24, .24, .04)}
    f, ax = plt.subplots(1, len(hm_matrices), figsize=(6, 5)) #, gridspec_kw=grid_kws)

    mask = np.zeros_like(hm_matrices[0], dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    for i in range(mask.shape[0]):
        mask[i,i] = False

    for ix, hm_matrix in enumerate(hm_matrices):

        hm_matrix = np.transpose(hm_matrix)
        
        if ix == 0:
            yticks = axis_labels
        else:
            yticks = False
        if ix == len(hm_matrices) - 1:
            cbar = False
            # cbar_ax = ax[len(hm_matrices)]
        else:
            cbar = False
            cbar_ax = None
        # f, ax = plt.subplots(figsize=(7, 5))
        ax[ix] = sns.heatmap(hm_matrix, cmap="Greens", mask=mask, vmin=0, center=0,
                ax=ax[ix],
                square=True, 
                linewidths=.3, 
                cbar=False,
                # cbar_ax=cbar_ax,
                # cbar_kws={"shrink": .1},
                annot=True, fmt='.2f', annot_kws={"size": 6},
                xticklabels=axis_labels, yticklabels=yticks)
        
        # ax[ix].set(title=titles[ix])
        ax[a,b].set_title(titles[ix], fontsize=8)

        plt.yticks(rotation=0)

    # plt.title(title)

    plt.tight_layout()
    # plt.show()
    if outfile:
        plt.savefig(outfile)

def is_same_cluster(classification, p1, p2):
    return classification[p1] == classification[p2]

def compile_pandas(simulations, tools, tool_names, measure_name):
    mat = np.zeros((simulations, len(tool_names)))

    for i in range(simulations):
        for j in range(len(tool_names)):
            mat[i][j] = tools[j][i]

    df = pd.DataFrame(mat, columns=tool_names).assign(Measure=measure_name)

    return df

def calc_prec_rec(ground, tool):
    precisions = list()
    recalls = list()
    for index, _ in enumerate(ground):
        clusters_ground = ground[index][0]
        clusters_tool = tool[index][0]

        tp = 0
        fp = 0
        fn = 0
        for key in clusters_ground:
            clust = clusters_ground[key]

            for pairs in itertools.combinations(clust, 2):
                if is_same_cluster(tool[index][1], *pairs):
                    tp += 1
                else:
                    fn += 1
        
        for key in clusters_tool:
            clust = clusters_tool[key]

            for pairs in itertools.combinations(clust, 2):
                if not is_same_cluster(ground[index][1], *pairs):
                    fp += 1
        
        prec = float(tp) / (tp + fp)
        rec = float(tp) / (tp + fn)

        precisions.append(prec)
        recalls.append(rec)

    return precisions, recalls

def plot_prec_and_rec(df_prec, df_rec, title=None, outfile=None):
    import seaborn as sns
    import matplotlib.pyplot as plt

    sns.set(style="whitegrid")

    f, ax = plt.subplots(1, 2, figsize=(9,3))

    ax[0] = sns.boxplot(data=df_prec, ax=ax[0], palette="Set2", linewidth=1.5, orient='h')
    # ax[0] = sns.swarmplot(data=df_prec, ax=ax[0], palette="Set2", alpha=.6, size=4, edgecolor="gray", linewidth=.5, orient='h')
    ax[0].set(title="Precision", xlim=(-0.05,1.05))

    ax[1] = sns.boxplot(data=df_rec, ax=ax[1], palette="Set2", linewidth=1.5, orient='h')
    # ax[1] = sns.swarmplot(data=df_rec, ax=ax[1], palette="Set2", alpha=.6, size=4, edgecolor="gray", linewidth=.5, orient='h')
    ax[1].set(title="Recall", xlim=(-0.05,1.05))
    ax[1].set_yticks([])

    if title:
        f.suptitle(title, fontsize=12)


    plt.tight_layout()
    
    if outfile:
        plt.savefig(outfile)

def plot_box_2x2(df_ari, df_fmi, df_vm, df_cs, title=None, outfile=None):
    import seaborn as sns
    import matplotlib.pyplot as plt

    sns.set(style="whitegrid")

    f, ax = plt.subplots(2, 2, figsize=(9,5))

    ax[0,0] = sns.boxplot(data=df_ari, ax=ax[0,0], palette="Set2", linewidth=1.5, orient='h')
    # ax[0] = sns.swarmplot(data=df_prec, ax=ax[0], palette="Set2", alpha=.6, size=4, edgecolor="gray", linewidth=.5, orient='h')
    ax[0,0].set(title='Adjusted Rand Index', xlim=(-0.05,1.05))

    ax[0,1] = sns.boxplot(data=df_fmi, ax=ax[0,1], palette="Set2", linewidth=1.5, orient='h')
    # ax[0] = sns.swarmplot(data=df_prec, ax=ax[0], palette="Set2", alpha=.6, size=4, edgecolor="gray", linewidth=.5, orient='h')
    ax[0,1].set(title='Fowlkes–Mallows Index', xlim=(-0.05,1.05))
    ax[0,1].set_yticks([])

    ax[1,0] = sns.boxplot(data=df_vm, ax=ax[1,0], palette="Set2", linewidth=1.5, orient='h')
    # ax[0] = sns.swarmplot(data=df_prec, ax=ax[0], palette="Set2", alpha=.6, size=4, edgecolor="gray", linewidth=.5, orient='h')
    ax[1,0].set(title='V-Measure Index', xlim=(-0.05,1.05))

    ax[1,1] = sns.boxplot(data=df_cs, ax=ax[1,1], palette="Set2", linewidth=1.5, orient='h')
    # ax[0] = sns.swarmplot(data=df_prec, ax=ax[0], palette="Set2", alpha=.6, size=4, edgecolor="gray", linewidth=.5, orient='h')
    ax[1,1].set(title='Completeness Score', xlim=(-0.05,1.05))
    ax[1,1].set_yticks([])

    if title:
        f.suptitle(title, fontsize=12)


    plt.tight_layout()
    
    if outfile:
        plt.savefig(outfile)

if __name__ == '__main__':
    
    try:
        exp = int(sys.argv[1])
    except:
        print('Usage: python3 clust_accuracies.py [EXPERIMENT_NUMBER]')
        sys.exit(1)


    folders = [
        '../data/exp%d/ground' % exp,
        '../data/exp%d/celluloid' % exp,
        '../data/exp%d/kmodes' % exp,
        '../data/exp%d/kmeans' % exp,
        '../data/exp%d/affinity' % exp,
        '../data/exp%d/agglomerative' % exp,
        '../data/exp%d/birch' % exp,
        '../data/exp%d/spectral' % exp
    ]
    filenames = [
        '%s/sim_%d_log_cluster_mutations.txt',
        '%s/sim_%d_scs_celluloid_clusters.txt',
        '%s/sim_%d_scs_kmodes_clusters.txt',
        '%s/sim_%d_scs_kmeans_clusters.txt',
        '%s/sim_%d_scs_affinity_clusters.txt',
        '%s/sim_%d_scs_agglomerative_clusters.txt',
        '%s/sim_%d_scs_birch_clusters.txt',
        '%s/sim_%d_scs_spectral_clusters.txt'
    ]
    cluster_labels = [
        'ground',
        'celluloid',
        'k-modes',
        'k-means',
        'affinity',
        'agglomerative',
        'BIRCH',
        'spectral'
    ]
    tot_simulations = 50

    import os
    if not os.path.exists('../plots/exp{0}'.format(exp)):
        os.makedirs('../plots/exp{0}'.format(exp))
    out_file = '../plots/exp{0}/exp{0}_plot_%s.pdf'.format(exp)

    assert(len(folders) == len(filenames) == len(cluster_labels))

    clusters = [[] for x in range(len(filenames))]

    for sim in range(1, tot_simulations+1):
        for c in range(len(clusters)):
            clusters[c].append(
                load_clusters.load_cluster_file(filenames[c] % (folders[c], sim))
            )

    # ---------------------------
    # --- Cluster Accuracies ----
    # ---------------------------

    mat_ari, fc_ari = build_heatmap_matrix(clusters, score=adjusted_rand_score)
    print('ARI')
    print(mat_ari)
    plot_heatmap(mat_ari, cluster_labels, 'Adjusted Rand Index', outfile=out_file % 'ari')

    mat_fmi, fc_fmi = build_heatmap_matrix(clusters, score=fowlkes_mallows_score)
    print('\nFMI')
    print(mat_fmi)
    plot_heatmap(mat_fmi, cluster_labels, 'Fowlkes–Mallows Index', outfile=out_file % 'fmi')

    mat_vm, fc_vm = build_heatmap_matrix(clusters, score=v_measure_score)
    print('\nVM')
    print(mat_vm)
    plot_heatmap(mat_vm, cluster_labels, 'V-Measure Index', outfile=out_file % 'vm')

    mat_cs, fc_cs = build_heatmap_matrix(clusters, score=completeness_score)
    print('\nCS')
    print(mat_cs)
    plot_heatmap(mat_cs, cluster_labels, 'Completeness Score', outfile=out_file % 'cs')

    plot_heatmaps_compacted_2x2([mat_ari, mat_fmi, mat_vm, mat_cs], cluster_labels, [
        'Adjusted Rand Index', 'Fowlkes–Mallows Index', 'V-Measure Index', 'Completeness Score'
    ], outfile=out_file % 'compact')

    # ---------------------------
    # ---- Precision-Recall -----
    # ---------------------------

    cluster_precs = list()
    cluster_recs = list()

    for c in range(len(clusters) - 1):
        prec, rec = calc_prec_rec(clusters[0], clusters[c+1])
        cluster_precs.append(prec)
        cluster_recs.append(rec)

    df_prec = compile_pandas(tot_simulations, 
                                cluster_precs, 
                                cluster_labels[1:],
                                'Precision')

    df_rec = compile_pandas(tot_simulations, 
                                cluster_recs, 
                                cluster_labels[1:],
                                'Recall')

    plot_prec_and_rec(df_prec, df_rec, outfile=out_file % 'pr')

    df_ari = compile_pandas(tot_simulations, 
                                fc_ari, 
                                cluster_labels[1:],
                                'Adjusted Rand Index')

    df_fmi = compile_pandas(tot_simulations, 
                                fc_fmi, 
                                cluster_labels[1:],
                                'Fowlkes–Mallows Index')

    df_vm = compile_pandas(tot_simulations, 
                                fc_vm, 
                                cluster_labels[1:],
                                'V-Measure Index')
    
    df_cs = compile_pandas(tot_simulations, 
                                fc_cs, 
                                cluster_labels[1:],
                                'Completeness Score')

    plot_box_2x2(df_ari, df_fmi, df_vm, df_cs, outfile=out_file % 'fcbox')