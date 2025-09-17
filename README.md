# ddmmseqs

## Description

<u>**D**</u>IAMOND-<u>**d**</u>ivide followed by <u>**mmseqs**</u>-linclust to cluster very large (>1B) protein sequence collections.

Uses a divide-and-conquer strategy to make clustering tractable. Relies on nature of the proteins being such that they will naturally split into distinct groups without any in-between sequences linking them, otherwise greedy clustering will just produce one big blob cluster.
* ensure that alignment coverage threshold is fairly high (>0.8) and this should be a valid assumption

## Pipeline overview

1. Build DIAMOND index from file-list (defined in samplesheet)
2. With DIAMOND bastp, query this protein database against itself to a secondary-AAI threshold
3. Pipe these results into a custom Python script that keeps track of greedy clusters
4. Gather sequences into a single file for each greedy-cluster
5. Chunk these greedy-clusters and cluster with `mmseqs linclust`

## Authors

Timothy J. Rozday