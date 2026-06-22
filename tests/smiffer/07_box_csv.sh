#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 7: Trajectory mode with per-frame box CSV"

folder_few_frames=testdata/smiffer/trajs/few_frames
path_pdb=$folder_few_frames/fse.pdb
path_traj=$folder_few_frames/fse.xtc
path_boxes=$folder_few_frames/boxes.csv

### per-frame box driven by a CSV (x_min x_max y_min y_max z_min z_max per row)
python3 volgrids smiffer $path_pdb -o $folder_few_frames -t $path_traj -B $path_boxes

### same, but with trajectory multiprocessing to exercise the forked-worker path
python3 volgrids smiffer $path_pdb -o $folder_few_frames -t $path_traj -B $path_boxes -n 4
