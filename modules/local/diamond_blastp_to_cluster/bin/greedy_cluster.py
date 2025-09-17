import sys
import os
from collections import defaultdict
import gzip
import itertools as it
import argparse

parser = argparse.ArgumentParser(description='Single linkage cluster based on pipe of PAF.')
parser.add_argument('-f', "--filelist", type=str,
                    help="Path to file containing list of .faa or .faa.gz files.")
parser.add_argument('-n', "--out_nodes", type=str,
                    help="Output nodes file")
parser.add_argument('-l', "--out_linkage", type=str,
                    help="Output linkage file")
parser.add_argument('-c', "--out_clusters", type=str,
                    help="Output clusters file.")
parser.add_argument('-s', "--out_cluster_seqs_dir", type=str,
                    help="Output directory of cluster faa files.")
parser.add_argument('-t', "--target_cluster_size", type=int, default=-1, help="Target number of members of a cluster for output files.")

args = parser.parse_args()


def parse_faa(lines):
    header = None
    seq_lines = []
    for l in lines:
        if l[0]=='>':
            if (header is not None) and seq_lines:
                yield (header, ''.join(seq_lines))
            header = l.strip()[1:].split()[0]
            seq_lines = []
        else:
            if l:
                seq_lines.append(l)
    else:
        if seq_lines:
            yield (header, ''.join(seq_lines))


cluster2nodes = defaultdict(set)
node2cluster = {}
node_index = {}
with gzip.open(args.out_nodes, 'wt') as out_nodes, gzip.open(args.out_linkage, 'wt') as out_linkage:
    for line in sys.stdin.readlines():
        line_ = [v.strip() for v in line.split()]
        node1, node2 = line_[:2]

        # add node to index and to clusters
        for node in [node1, node2]:
            if not node in node_index:
                node_index[node] = len(node_index)
                out_nodes.write(f"{node}\t{node_index[node]}\n")

            if not node in node2cluster:
                cluster_idx = node_index[node]
                node2cluster[node_index[node]] = cluster_idx
                cluster2nodes[cluster_idx] = {node_index[node]}
        
        # add linkage to file
        out_linkage.write(f"{node_index[node1]}\t{node_index[node2]}\n")

        # create/merge clusters
        node1_idx, node2_idx = node_index[node1], node_index[node2] 
        if node2cluster[node1_idx] == node2cluster[node2_idx]:
            continue

        keep_cluster_idx = node1_idx
        delete_cluster_idx = node2_idx
        if node1_idx > node2_idx:
            keep_cluster_idx, delete_cluster_idx = delete_cluster_idx, keep_cluster_idx
        
        cluster2nodes[keep_cluster_idx].update(cluster2nodes[delete_cluster_idx])
        
        for n in cluster2nodes[delete_cluster_idx]:
            node2cluster[n] = keep_cluster_idx
        del cluster2nodes[delete_cluster_idx]


# write cluster information
with gzip.open(args.out_clusters, 'wt') as f:
    for k,v in node2cluster.items():
        f.write(f"{k}\t{v}\n")


# write cluster sequences

# first merge clusters up to target size
if args.target_cluster_size>0:
    group2clusters = defaultdict(set)
    cluster2group = {}

    cluster_sizes = {k:len(vs) for k,vs in cluster2nodes.items()}
    sorted_cluster_sizes = list(sorted(cluster_sizes.items(), key=lambda x:-x[1]))
    for k1,_ in sorted_cluster_sizes:
        if k1 in cluster2group:
            continue

        cluster2group[k1] = len(group2clusters)
        group2clusters[cluster2group[k1]] = {k1}

        group_size = sum([cluster_sizes[v] for v in group2clusters[cluster2group[k1]]])
        if group_size>=args.target_cluster_size:
            continue

        for k2,s in sorted_cluster_sizes:
            if k2 in cluster2group:
                continue

            if group_size+s <= args.target_cluster_size:
                cluster2group[k2] = cluster2group[k1]
                group2clusters[cluster2group[k1]].add(k2)
                group_size += s
else:
    group2clusters = {i:{k} for i,k in enumerate(cluster2nodes)}
    cluster2group = {v:k for k,vs in group2clusters.items() for v in vs}


node2group = {k:cluster2group[v] for k,v in node2cluster.items()}
group2nodes = defaultdict(set)
for k,v in node2group.items():
    group2nodes[v].add(k)

os.makedirs(args.out_cluster_seqs_dir, exist_ok=True)
cluster_files = {k: gzip.open(f"{args.out_cluster_seqs_dir}/{k}.faa.gz", 'wt')
                 for k,_ in group2nodes.items()}

with open(args.filelist, 'rt') as f_list:
    for l in f_list:
        fp = l.strip()
        file_opener = gzip.open if fp.endswith('.gz') else open

        with file_opener(fp, 'rt') as f_faa:
            for header, seq in parse_faa(l.strip() for l in f_faa.readlines()):
                node_idx = node_index[header]
                cluster_idx = node2group[node_idx]
                cluster_files[cluster_idx].write(f">{header}\n{seq}\n\n")

for k,f in cluster_files.items():
    f.close()


            
