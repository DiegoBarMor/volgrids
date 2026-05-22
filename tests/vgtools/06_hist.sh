#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 6: Voxel histribution histogram (vgtools histogram)"

fhist="testdata/vgtools/hist"
fc="testdata/vgtools/converting"
ff="testdata/vgtools/fix_cmap"
fu="testdata/vgtools/unpacking"

mkdir -p "$fhist"


############################# SINGLE-FRAME DX
echo ">>> hist: single-frame DX"
python3 volgrids vgtools histogram "$fc/1iqj.stk.dx" --out "$fhist/hist_dx.png"


############################# SINGLE-FRAME CMAP
echo ">>> hist: single-frame CMAP"
python3 volgrids vgtools histogram "$fc/1iqj.stk.cmap" --out "$fhist/hist_cmap_single.png"


############################# MULTI-FRAME CMAP (all frames)
echo ">>> hist: multi-frame CMAP (all frames)"
python3 volgrids vgtools histogram "$ff/hbdonors.issue.cmap" --out "$fhist/hist_cmap_all.png"


############################# MULTI-FRAME CMAP (single key)
echo ">>> hist: multi-frame CMAP (single key)"
python3 volgrids vgtools histogram "$fu/1iqj.cmap" --key "1iqj.hbacceptors" --out "$fhist/hist_cmap_key.png"

# [TODO] add test for grid with all 0 values
