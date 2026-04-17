#!/bin/bash
set -eu

folder_env="testdata/env"
folder_smiffer="testdata/smiffer"
folder_vgtools="testdata/vgtools"

rm -rf $folder_env

# rm -rf "$folder_smiffer/apbs"
rm  -f "$folder_smiffer/toy_systems/"*.cmap
rm -rf "$folder_smiffer/pocket_sphere"
rm -rf "$folder_smiffer/whole"*
rm  -f "$folder_smiffer/trajs/few_frames/"*.cmap
rm  -f "$folder_smiffer/trajs/few_resids/"*.cmap
rm  -f "$folder_smiffer/ligands/all_interactions/"*.cmap
rm  -f "$folder_smiffer/ligands/mid_size/"*.cmap
rm  -f "$folder_smiffer/ligands/nr3_case/"*.cmap
rm -rf "$folder_smiffer/cavities"

folder_vgt_00="$folder_vgtools/converting"
folder_vgt_01="$folder_vgtools/packing"
folder_vgt_02="$folder_vgtools/unpacking"
folder_vgt_03="$folder_vgtools/fix_cmap"
rm -f $folder_vgt_00/dx* $folder_vgt_00/mrc* $folder_vgt_00/ccp4* $folder_vgt_00/cmap*
rm -f $folder_vgt_01/2esj.cmap
rm -f $folder_vgt_02/1iqj.*.cmap
rm -f $folder_vgt_03/hbdonors.fixed.cmap
