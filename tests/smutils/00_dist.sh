#!/bin/bash
set -eu

echo
echo ">>> TEST SMUTILS 0: Voxel distribution (smutils dist)"

fdist="testdata/smutils/dist"
fc="testdata/vgtools/converting"
ff="testdata/vgtools/fix_cmap"
fu="testdata/vgtools/unpacking"

mkdir -p "$fdist"


############################# SINGLE-FRAME DX
echo ">>> dist: single-frame DX"
python3 volgrids smutils plot_dist "$fc/1iqj.stk.dx" --out "$fdist/dist_dx.png"


############################# SINGLE-FRAME CMAP
echo ">>> dist: single-frame CMAP"
python3 volgrids smutils plot_dist "$fc/1iqj.stk.cmap" --out "$fdist/dist_cmap_single.png"


############################# MULTI-FRAME CMAP — all frames
echo ">>> dist: multi-frame CMAP (all frames)"
python3 volgrids smutils plot_dist "$ff/hbdonors.issue.cmap" --out "$fdist/dist_cmap_all.png"


############################# MULTI-FRAME CMAP — single key
echo ">>> dist: multi-frame CMAP (single key)"
python3 volgrids smutils plot_dist "$fu/1iqj.cmap" --key "1iqj.hbacceptors" --out "$fdist/dist_cmap_key.png"
