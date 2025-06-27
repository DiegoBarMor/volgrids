#!/bin/bash
set -eu

folder_smiffer="testdata/smiffer"
folder_vgtools="testdata/vgtools"

folder00="$folder_smiffer/pocket_sphere"
folder01="$folder_smiffer/whole"
folder02="$folder_smiffer/traj"
folder03c="$folder_vgtools/converting"
folder03p="$folder_vgtools/packing"
folder03u="$folder_vgtools/unpacking"
folder03f="$folder_vgtools/fix_cmap"
folder04="$folder_smiffer/ligand"

rm -rf $folder00 $folder01
rm  -f $folder02/*.cmap $folder02/*.json

rm -f $folder03c/dx* $folder03c/mrc* $folder03c/cmap*
rm -f $folder03p/2esj.cmap
rm -f $folder03u/1iqj.*.cmap
rm -f $folder03f/hbdonors.fixed.cmap

rm -f $folder04/*.cmap $folder04/*.json

clear
