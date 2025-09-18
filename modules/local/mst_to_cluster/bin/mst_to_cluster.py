import os
from collections import defaultdict
import gzip
import datetime
import argparse

parser = argparse.ArgumentParser(description='Single linkage cluster based on pipe of PAF.')
parser.add_argument('-l', "--linkagelists", type=str,
                    help="Path to TSV file containing paths to node indexes and linkages.")
parser.add_argument('-f', "--filelists", type=str,
                    help="Path to file containing list of paths to files containing lists of .faa or .faa.gz files.")
parser.add_argument('-n', "--out_nodes", type=str,
                    help="Output nodes file")
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
master_node_index = {}
next_cluster_idx = 0
next_node_idx = 0
line_i = 0
with open(args.linkagelists, 'rt') as linkagelist_f, gzip.open(args.out_nodes, 'wt') as out_nodes:
    for linkagelist_l in linkagelist_f:
        node_fp, linkage_fp = [v.strip() for v in linkagelist_l.split()]

        node_index = {}
        with gzip.open(node_fp, 'rt') as node_f:
            for node_l in node_f:
                k, v = [v.strip() for v in node_l.split()]
                node_index[k] = int(v)
        rev_node_index = {}
        for k,v in node_index.items():
            rev_node_index[v] = k
            if not k in master_node_index:
                master_node_index[k] = next_node_idx
                out_nodes.write(f"{k}\t{master_node_index[k]}\n")
                next_node_idx += 1
            if not master_node_index[k] in node2cluster:
                cluster_idx = next_cluster_idx
                node2cluster[master_node_index[k]] = cluster_idx
                cluster2nodes[cluster_idx] = {master_node_index[k]}
                next_cluster_idx += 1

        with gzip.open(linkage_fp, 'rt') as linkage_f:
            for linkage_l in linkage_f:
                line_i += 1
                if line_i % 1_000_000:
                    print(f"{datetime.datetime.now()}\t{line_i:,} lines read")
                node1, node2 = [rev_node_index[int(v.strip())] for v in linkage_l.split()]
                node1_index, node2_index = master_node_index[node1], master_node_index[node2]
                
                # create/merge clusters
                keep_cluster_idx, delete_cluster_idx = node2cluster[node1_index], node2cluster[node2_index]
                if keep_cluster_idx == delete_cluster_idx:
                    continue

                if keep_cluster_idx > delete_cluster_idx:
                    keep_cluster_idx, delete_cluster_idx = delete_cluster_idx, keep_cluster_idx
        
                cluster2nodes[keep_cluster_idx].update(cluster2nodes[delete_cluster_idx])
        
                for n in cluster2nodes[delete_cluster_idx]:
                    node2cluster[n] = keep_cluster_idx
                del cluster2nodes[delete_cluster_idx]


# write cluster information
with gzip.open(args.out_clusters, 'wt') as f:
    for k,v in node2cluster.items():
        f.write(f"{k}\t{v}\n")


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


# write cluster sequences
node2group = {k:cluster2group[v] for k,v in node2cluster.items()}
group2nodes = defaultdict(set)
for k,v in node2group.items():
    group2nodes[v].add(k)

os.makedirs(args.out_cluster_seqs_dir, exist_ok=True)
cluster_files = {k: gzip.open(f"{args.out_cluster_seqs_dir}/{k}.faa.gz", 'wt')
                 for k,_ in group2nodes.items()}

seq_i = 0
with open(args.filelists, 'rt') as f_lists:
    for filelists_l in f_lists:
        filelist_fp = filelists_l.strip()
        with open(filelist_fp, 'rt') as f_list:
            for l in f_list:
                fp = l.strip()
                file_opener = gzip.open if fp.endswith('.gz') else open

                with file_opener(fp, 'rt') as f_faa:
                    for header, seq in parse_faa(l.strip() for l in f_faa.readlines()):
                        seq_i += 1
                        if seq_i % 1_000_000:
                            print(f"{datetime.datetime.now()}\t{seq_i:,} seqs written")
                        node_idx = master_node_index[header]
                        cluster_idx = node2group[node_idx]
                        cluster_files[cluster_idx].write(f">{header}\n{seq}\n\n")

for k,f in cluster_files.items():
    f.close()


            
