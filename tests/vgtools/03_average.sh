#!/bin/bash
set -euo pipefail

echo
echo ">>> TEST VGTOOLS 3: Average of trajectory grids"

folder_in="testdata/smiffer/trajs/small_rna"
folder_out="testdata/vgtools/average"
conf_just_stacking="SMIF_STK=True SMIF_HBA=False SMIF_HBD=False SMIF_HPHOB=False SMIF_HPHIL=False SMIF_APBS=False TRIM_SAVE=False"

mkdir -p "$folder_out"

python3 volgrids smiffer "$folder_in/struct_A.pdb" -o "$folder_out" -t "$folder_in/minA.xtc" --config "$conf_just_stacking"

python3 volgrids vgtools average "$folder_out/struct_A.stk.cmap" "$folder_out/struct_A.stk.avg.cmap"
