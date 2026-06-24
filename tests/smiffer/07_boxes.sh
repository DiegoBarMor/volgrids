#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 7: Trajectory mode with custom boxes"

folder_few_frames=testdata/smiffer/trajs/few_frames
path_pdb=$folder_few_frames/fse.pdb
path_traj=$folder_few_frames/fse.xtc
path_boxes=$folder_few_frames/boxes.csv

conf_just_stak="GRID_FORMAT_OUTPUT=CMAP DO_SMIF_STACKING=True DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False"

### per-frame box driven by a CSV (x_min x_max y_min y_max z_min z_max per row)
python3 volgrids smiffer $path_pdb -o $folder_few_frames -t $path_traj -c "$conf_just_stak" -B $path_boxes

### same, but with trajectory multiprocessing to exercise the forked-worker path
python3 volgrids smiffer $path_pdb -o $folder_few_frames -t $path_traj -c "$conf_just_stak" -B $path_boxes -n 4
