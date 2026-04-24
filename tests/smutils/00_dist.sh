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


############################# 3D VOLUME — single-frame DX
echo ">>> plot_3d: single-frame DX"
python3 volgrids smutils plot_3d "$fc/1iqj.stk.dx" --out "$fdist/vol_dx.html" 


############################# 3D VOLUME — multi-frame CMAP averaged
echo ">>> plot_3d: multi-frame CMAP (averaged)"
python3 volgrids smutils plot_3d "$ff/hbdonors.issue.cmap" --out "$fdist/vol_cmap_avg.html" 


############################# 3D VOLUME — single key
echo ">>> plot_3d: multi-frame CMAP (single key)"
python3 volgrids smutils plot_3d "$fu/1iqj.cmap" --key "1iqj.hbacceptors" --out "$fdist/vol_cmap_key.html"


############################# ECHARTS — single-frame DX
echo ">>> plot_echarts: single-frame DX"
python3 volgrids smutils plot_echarts "$fc/1iqj.stk.dx" --out "$fdist/echarts_dx.html"


############################# ECHARTS — multi-frame CMAP averaged
echo ">>> plot_echarts: multi-frame CMAP (averaged)"
python3 volgrids smutils plot_echarts "$ff/hbdonors.issue.cmap" --out "$fdist/echarts_cmap_avg.html"


############################# ECHARTS — single key
echo ">>> plot_echarts: multi-frame CMAP (single key)"
python3 volgrids smutils plot_echarts "$fu/1iqj.cmap" --key "1iqj.hbacceptors" --out "$fdist/echarts_cmap_key.html"
