#!/bin/bash
set -euo pipefail

fdata="testdata/smiffer/trajs/spheres"
path_pdb=$fdata/rna.pdb
path_traj=$fdata/rna.xtc
radius_extra=0.0 # override as needed (default is 0.0)

calc_spheres() {
    idx=$1
    query=$2
    fout="$fdata/sphere_$idx"
    echo "Calculating spheres for query: $query"
    info_spheres=$(python3 volgrids smutils sphere find $path_pdb "$query" -t $path_traj -r $radius_extra)
    python3 volgrids smutils sphere grid $path_pdb -s "$info_spheres" -t $path_traj -o "$fout"
}
calc_spheres 1 "resid 7 8" &\
calc_spheres 2 "resid 9 10" &\
calc_spheres 3 "resid 11 12" &\
calc_spheres 4 "resid 20 21" &\
calc_spheres 5 "resid 22 23" &\
calc_spheres 6 "resid 24 25"
wait


merge() {
    python3 volgrids vgtools op or "$fdata/$1/rna.sphere.cmap" "$fdata/$2/rna.sphere.cmap" "$fdata/$3"
}
merge "sphere_1" "sphere_2" "merged_1_2.cmap" &\
merge "sphere_3" "sphere_4" "merged_3_4.cmap" &\
merge "sphere_5" "sphere_6" "merged_5_6.cmap"
wait
rm -rf "$fdata/sphere_"*


merge() {
    python3 volgrids vgtools op or "$fdata/$1" "$fdata/$2" "$fdata/$3"
}
merge "merged_1_2.cmap" "merged_3_4.cmap" "merged_1_4.cmap"
merge "merged_1_4.cmap" "merged_5_6.cmap" "rna.sphere.cmap"
rm -rf "$fdata/merged_"*
