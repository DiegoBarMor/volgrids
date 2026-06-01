#!/bin/bash
set -euo pipefail

##### test with simple toy system SMIF
# fout="testdata/smiffer/toy_systems"
# conf="GRID_FORMAT_OUTPUT=BIN DO_SMIF_STACKING=True DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False"
# python3 volgrids smiffer $fout/guanine.pdb -o $fout -c $conf

# path_smif="$fout/guanine.stacking.bin"
# path_clusters="$fout/guanine.stacking.clusters.bin"
# threshold=0.01


# python3 volgrids vgtools convert "$path_smif" -m


# folder_src="examples/bin_format"
# gcc -O2 $folder_src/clusters.c -o $folder_src/clusters
# ./$folder_src/clusters "$path_smif" "$path_clusters" $threshold
# python3 volgrids vgtools convert "$path_clusters" -m


##### test with larger SMIF (must run examples/interfaces/run.sh first)
fdata="testdata/smiffer/interfaces/prot_rna"
ftmp="$fdata/tmp_unpack"

path_smif="$fdata/prot.stacking.bin"
path_clusters="$fdata/prot.stacking.clusters.bin"
threshold=0.1

mkdir -p "$ftmp"
python3 volgrids vgtools unpack "$fdata/prot.smif.cmap" "$ftmp"
python3 volgrids vgtools convert "$ftmp/prot.stacking.cmap" -b "$path_smif"
rm -rf "$ftmp"

folder_src="examples/bin_format"
gcc -O2 $folder_src/clusters.c -o $folder_src/clusters
./$folder_src/clusters "$path_smif" "$path_clusters" $threshold

python3 volgrids vgtools convert "$path_clusters" -m
