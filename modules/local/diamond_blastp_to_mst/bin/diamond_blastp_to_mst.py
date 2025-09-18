import sys
import os
from collections import defaultdict
import gzip
import datetime
import argparse

parser = argparse.ArgumentParser(description='Single linkage cluster based on pipe of PAF.')
parser.add_argument('-n', "--out_nodes", type=str,
                    help="Output nodes file")
parser.add_argument('-l', "--out_linkage", type=str,
                    help="Output linkage file")
args = parser.parse_args()


cluster2nodes = defaultdict(set)
node2cluster = {}
node_index = {}
next_cluster_idx = 0
next_node_idx = 0
line_i = 0
with gzip.open(args.out_nodes, 'wt') as out_nodes, gzip.open(args.out_linkage, 'wt') as out_linkage:
    for line in sys.stdin.readlines():
        line_i += 1
        if line_i % 1_000_000:
            print(f"{datetime.datetime.now()}\t{line_i:,} lines read")

        line_ = [v.strip() for v in line.split()]

        # self-match lines are truncated        
        if len(line_)==11:
            node  = line_[0]

            if not node in node_index:
                node_index[node] = next_node_idx
                next_node_idx += 1
                out_nodes.write(f"{node}\t{node_index[node]}\n")

            if not node_index[node] in node2cluster:
                cluster_idx = next_cluster_idx
                node2cluster[node_index[node]] = cluster_idx
                cluster2nodes[cluster_idx] = {node_index[node]}
                next_cluster_idx += 1

            continue

        query_name, query_length, query_start, \
        query_end, strand, target_name, \
        target_length, target_start, target_end, \
        n_matching, n_bases, mapping_quality, \
        bit_score, raw_score, evalue = line_

        node1, node2 = query_name, target_name

        # add node to index and to clusters
        for node in [node1, node2]:
            if not node in node_index:
                node_index[node] = next_node_idx
                next_node_idx += 1
                out_nodes.write(f"{node}\t{node_index[node]}\n")

            if not node_index[node] in node2cluster:
                cluster_idx = next_cluster_idx
                node2cluster[node_index[node]] = cluster_idx
                cluster2nodes[cluster_idx] = {node_index[node]}
                next_cluster_idx += 1
        
        # create/merge clusters
        keep_cluster_idx, delete_cluster_idx = node2cluster[node_index[node1]], node2cluster[node_index[node2]]
        if keep_cluster_idx == delete_cluster_idx:
            continue

        # add linkage to file
        out_linkage.write(f"{node_index[node1]}\t{node_index[node2]}\n")
        
        if keep_cluster_idx > delete_cluster_idx:
            keep_cluster_idx, delete_cluster_idx = delete_cluster_idx, keep_cluster_idx
        
        cluster2nodes[keep_cluster_idx].update(cluster2nodes[delete_cluster_idx])
        
        for n in cluster2nodes[delete_cluster_idx]:
            node2cluster[n] = keep_cluster_idx
        del cluster2nodes[delete_cluster_idx]
