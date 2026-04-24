#!/bin/bash
set -eu

folder_env="testdata/env"
folder_smiffer="testdata/smiffer"
folder_smutils="testdata/smutils"
folder_vgtools="testdata/vgtools"

rm -rf $folder_env

# rm -rf "$folder_smiffer/apbs"
rm  -f "$folder_smiffer/toy_systems/"*.cmap
rm -rf "$folder_smiffer/pocket_sphere"
rm -rf "$folder_smiffer/whole"*
rm  -f "$folder_smiffer/trajs/"**/*.cmap
rm  -f "$folder_smiffer/ligands/"**/*.cmap
rm  -f "$folder_smiffer/interfaces/"**/*.cmap
rm -rf "$folder_smiffer/cavities"

rm -rf "$folder_smutils/hist"

rm -f $folder_vgtools/converting/dx*
rm -f $folder_vgtools/converting/mrc*
rm -f $folder_vgtools/converting/ccp4*
rm -f $folder_vgtools/converting/cmap*
rm -f $folder_vgtools/packing/2esj.cmap
rm -f $folder_vgtools/unpacking/1iqj.*.cmap
rm -f $folder_vgtools/fix_cmap/hbdonors.fixed.cmap
rm -f $folder_vgtools/operations/*.cmap
