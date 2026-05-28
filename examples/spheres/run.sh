#!/bin/bash
set -euo pipefail

fdata="testdata/smiffer/trajs/spheres"
path_pdb=$fdata/struct.pdb
path_traj=$fdata/traj.xtc
radius_extra=2.0 # override as needed (default is 2.0)

calc_spheres() {
    query=$1
    echo "Calculating spheres for query: $query"
    info_spheres=$(python3 volgrids smutils sphere find $path_pdb "$query" -t $path_traj -r $radius_extra)
    python3 volgrids smutils sphere grid $path_pdb -s $info_spheres -t $path_traj
}

calc_spheres "resid 7 8"
