#!/bin/bash
set -euo pipefail

path_smif=$1
path_clusters=$2
threshold=$3

folder_src="examples/bin_format"
gcc -O2 $folder_src/clusters.c -o $folder_src/clusters
./$folder_src/clusters "$path_smif" "$path_clusters" "$threshold"

python3 examples/bin_format/sort_clusters.py "$path_clusters" 200
python3 examples/bin_format/expand_mrc2cmap.py "$path_clusters" # outputs a CMAP

rm -f "$path_clusters"
