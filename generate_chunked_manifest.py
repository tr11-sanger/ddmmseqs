import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Chunk input filelist and generate manifest for ddmmseqs Nextflow pipeline')
parser.add_argument('-f', '--filelist', type=str, help='Path to file containing list of FAA files.')
parser.add_argument('-c', '--name', type=str, help='Name of protein sequence collection.')
parser.add_argument('-n', '--num_chunks', type=int, help='Number of chunks.')
parser.add_argument('-p', '--chunk_prefix', type=str, help='Prefix for output filelist chunks.')
parser.add_argument('-o', '--output', type=str, help='Path to output manifest CSV file.')
args = parser.parse_args()

with open(args.filelist, 'rt') as f_in:
    fps = [v.strip() for v in f_in.readlines()]

out_fps = {i:f"{args.chunk_prefix}_{i}.txt" for i in range(1,args.num_chunks+1)}

file_write_buffer = defaultdict(list)
for i,fp in enumerate(fps):
    out_fp = out_fps[int(((i/len(fps))*args.num_chunks))+1]
    file_write_buffer[out_fp].append(f"{fp}\n")

for out_fp, lines in file_write_buffer.items():
    with open(out_fp, 'wt') as f_out:
        f_out.writelines(lines)


manifest = []
for chunk_n, out_fp in out_fps.items():
    manifest.append({
        'collection_name': args.name,
        'chunk_name': f"chunk_{chunk_n}",
        'filelist': out_fp,
    })
manifest = pd.DataFrame(manifest)
manifest.to_csv(args.output, index=False)
