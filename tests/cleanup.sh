#!/bin/bash
set -eu

folder_env="testdata/env"
folder_smiffer="testdata/smiffer"
folder_vgtools="testdata/vgtools"

folder00="$folder_smiffer/toy_systems"
folder01="$folder_smiffer/pocket_sphere"
folder02="$folder_smiffer/whole"
folder03="$folder_smiffer/trajs"
folder04="$folder_smiffer/ligands"
folder05c="$folder_vgtools/converting"
folder05p="$folder_vgtools/packing"
folder05u="$folder_vgtools/unpacking"
folder05f="$folder_vgtools/fix_cmap"

rm -rf $folder_env $folder00 $folder01 $folder02 $folder03 $folder04
rm -f $folder05c/dx* $folder05c/mrc* $folder05c/ccp4* $folder05c/cmap*
rm -f $folder05p/2esj.cmap
rm -f $folder05u/1iqj.*.cmap
rm -f $folder05f/hbdonors.fixed.cmap
