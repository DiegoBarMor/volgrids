#!/bin/bash
set -euo pipefail

##### test with simple toy system SMIF
fout="testdata/smiffer/toy_systems"
conf="OUT_FORMAT=BIN SMIF_STK=True SMIF_HBA=False SMIF_HBD=False SMIF_HPHOB=False SMIF_HPHIL=False SMIF_APBS=False"
python3 volgrids smiffer $fout/guanine.pdb -o $fout -c "$conf"

path_smif="$fout/guanine.stk.bin"
path_clusters="$fout/guanine.stk.clusters.bin"
isovalue=0.01

python3 volgrids vgtools convert "$path_smif" -f MRC
python3 volgrids vgtools segment "$path_smif" "$path_clusters" -i $isovalue

python3 examples/bin_format/expand_mrc2cmap.py "$path_clusters" # outputs a CMAP
rm -f "$path_clusters"
