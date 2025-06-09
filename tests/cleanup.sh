#!/bin/bash
set -eu

folder_data="data/tests"
folder00="$folder_data/00-pocket_sphere"
folder01="$folder_data/01-whole"
folder02="$folder_data/02-traj"
folder03="$folder_data/03-tools"

folder03c="$folder03/converting"
folder03p="$folder03/packing"
folder03u="$folder03/unpacking"
folder03f="$folder03/fix_cmap"

rm -rf $folder00 $folder01
rm  -f $folder02/*.cmap $folder02/*.json

rm -f $folder03c/dx* $folder03c/mrc* $folder03c/cmap*
rm -f $folder03p/2esj.cmap
rm -f $folder03u/1iqj.*.cmap
rm -f $folder03f/hbdonors.fixed.cmap
