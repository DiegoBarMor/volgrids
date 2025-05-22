#!/bin/bash
set -eu

folder_data="data/tests"
folder00="$folder_data/00-ps"
folder01="$folder_data/01-whole"
folder02="$folder_data/02-conversions"
folder03="$folder_data/03-cmap"
folder04="$folder_data/04-traj"

rm -rf $folder00 $folder01
rm -f $folder02/dx* $folder02/mrc*
rm -f $folder03/1iqj.*.mrc $folder03/2esj.cmap
