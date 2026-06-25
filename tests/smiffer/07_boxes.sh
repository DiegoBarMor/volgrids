#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 7: Trajectory mode with custom boxes"

folder_few_frames=testdata/smiffer/trajs/few_frames
path_pdb=$folder_few_frames/fse.pdb
path_traj=$folder_few_frames/fse.xtc
path_csv=$folder_few_frames/boxes.csv
path_txt=$folder_few_frames/boxes.txt

conf_just_stk="GRID_FORMAT_OUTPUT=CMAP DO_SMIF_STACKING=True DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False"


###### Priority 1: Box to be used is specifically requested by the user (via CLI)
### per-frame box driven by a CSV (x_min x_max y_min y_max z_min z_max per row)
python3 volgrids smiffer $path_pdb -t $path_traj -o $folder_few_frames/01_csv -c "$conf_just_stk" -b $path_csv

### same, but with trajectory multiprocessing to exercise the forked-worker path
python3 volgrids smiffer $path_pdb -t $path_traj -n 4 -o $folder_few_frames/01_csv_mp -c "$conf_just_stk" -b $path_csv

### box information provided as space-separated numbers fed directly to the -b flag
python3 volgrids smiffer $path_pdb -t $path_traj -n 4 -o $folder_few_frames/01_cli -c "$conf_just_stk" -b "$(cat $path_txt)"


###### Priority 2: Box to be used is based on a user-provided sphere (via CLI)
query="resid 35 36 37 38"
spheres=$(python3 volgrids smutils  sphere find $path_pdb "$query" -t $path_traj -r 2.0)
python3 volgrids smiffer $path_pdb -t $path_traj -n 4 -o $folder_few_frames/02_cli -c "$conf_just_stk" -s "$spheres"


##### Priority 3: dealing with a trajectory and the user requested the box to be inferred from the structure at every frame (via config `BOX_TIGHT_TRAJ=True`)
python3 volgrids smiffer $path_pdb -t $path_traj -n 4 -o $folder_few_frames/03_tight -c "$conf_just_stk BOX_TIGHT_TRAJ=True" -b $path_csv


##### Priority 4 (default): a common box that can enclose the structure in all frames is used (via config `BOX_TIGHT_TRAJ=False`)
python3 volgrids smiffer $path_pdb -t $path_traj -n 4 -o $folder_few_frames/04_common -c "$conf_just_stk BOX_TIGHT_TRAJ=False" -b $path_csv
