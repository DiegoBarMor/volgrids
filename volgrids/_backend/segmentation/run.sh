#!/bin/bash
set -euo pipefail

### This backend script is meant to be used by the Python scripts.
### It's not intended to be called directly by the user
### If you choose to run it directly, make sure to check the arguments are appropriate.

folder_src=$(realpath "$(dirname "$0")")

if [[ "$#" -ne 4 ]]; then
    echo "Usage: $0 <path_grid> <path_clusters> <isovalue> <volume_threshold>"
    exit 1
fi

path_grid=$1
path_clusters=$2
isovalue=$3
volume_threshold=$4

if [[ ! -f "$folder_src/main" ]]; then
    "$folder_src/compile.sh"
fi

"$folder_src/main" "$path_grid" "$path_clusters" "$isovalue"

python3 "$folder_src/sort_clusters.py" "$path_clusters" "$volume_threshold"
