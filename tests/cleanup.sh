#!/bin/bash
set -eu

folder_smiffer="testdata/smiffer"
folder_vgtools="testdata/vgtools"

folder00="$folder_smiffer/toy_systems"
folder01="$folder_smiffer/pocket_sphere"
folder02="$folder_smiffer/whole"
folder03="$folder_smiffer/traj"
folder04c="$folder_vgtools/converting"
folder04p="$folder_vgtools/packing"
folder04u="$folder_vgtools/unpacking"
folder04f="$folder_vgtools/fix_cmap"
folder05="$folder_smiffer/ligand"

rm -rf $folder00 $folder01 $folder02
rm  -f $folder03/*.cmap

rm -f $folder04c/dx* $folder04c/mrc* $folder04c/ccp4* $folder04c/cmap*
rm -f $folder04p/2esj.cmap
rm -f $folder04u/1iqj.*.cmap
rm -f $folder04f/hbdonors.fixed.cmap

rm -f $folder05/*.cmap

clear
