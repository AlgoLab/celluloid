def load_cluster_file(file_path):
    # from collections import defaultdict

    # clusters = defaultdict(list)
    clusters = dict()
    point_to_clust = dict()
    with open(file_path, 'r') as fin:
        for line in fin:
            clust, points = line.strip().split()
            
            clust = int(clust)
            points = [int(x) for x in points.replace('"', '').split(',')]
            
            clusters[clust] = points

            for p in points:
                point_to_clust[p] = clust

    point_prediction = [point_to_clust[x] for x in sorted(point_to_clust)]

    return clusters, point_to_clust, point_prediction

if __name__ == '__main__':
    load_cluster_file('../test/sim_1_scs_kmean_clusters_cells.txt')