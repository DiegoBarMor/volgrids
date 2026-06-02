#!/bin/bash
set -euo pipefail

##### test with simple toy system SMIF
fout="testdata/smiffer/toy_systems"
conf="GRID_FORMAT_OUTPUT=BIN DO_SMIF_STACKING=True DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False"
python3 volgrids smiffer $fout/guanine.pdb -o $fout -c "$conf"

path_smif="$fout/guanine.stacking.bin"
path_clusters="$fout/guanine.stacking.clusters.bin"
threshold=0.01

python3 volgrids vgtools convert "$path_smif" -m

bash examples/bin_format/segmentation.sh "$path_smif" "$path_clusters" $threshold
