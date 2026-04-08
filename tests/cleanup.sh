#!/bin/bash
set -eu

folder_env="testdata/env"
folder_smiffer="testdata/smiffer"
folder_vgtools="testdata/vgtools"

folder_sm_00="$folder_smiffer/toy_systems"
folder_sm_01="$folder_smiffer/pocket_sphere"
folder_sm_02="$folder_smiffer/whole"
folder_sm_03="$folder_smiffer/trajs"
folder_sm_04="$folder_smiffer/ligands"
folder_sm_05="$folder_smiffer/cavities"

folder_vgt_00="$folder_vgtools/converting"
folder_vgt_01="$folder_vgtools/packing"
folder_vgt_02="$folder_vgtools/unpacking"
folder_vgt_03="$folder_vgtools/fix_cmap"

rm -rf $folder_env $folder_sm_00 $folder_sm_01 $folder_sm_02 $folder_sm_03 $folder_sm_04 $folder_sm_05
rm -f $folder_vgt_00/dx* $folder_vgt_00/mrc* $folder_vgt_00/ccp4* $folder_vgt_00/cmap*
rm -f $folder_vgt_01/2esj.cmap
rm -f $folder_vgt_02/1iqj.*.cmap
rm -f $folder_vgt_03/hbdonors.fixed.cmap
