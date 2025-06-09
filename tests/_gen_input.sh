#!/bin/bash
set -eu

fpdb="data/pdb"
ftests="data/tests/03-tools"
fc="$ftests/converting"
fp="$ftests/packing"
fu="$ftests/unpacking"
ff="$ftests/fix_cmap"

mkdir -p $fc $fp $fu $ff
rm -f $fc/* $fp/* $fu/*  $ff/*.cmap

args_formats="DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False DO_OUTPUT_DX=True DO_OUTPUT_MRC=True DO_OUTPUT_CMAP=True"
args_smifs="DO_SMIF_HBA=True DO_SMIF_HBD=True DO_SMIF_HYDROPHOBIC=True DO_SMIF_HYDROPHILIC=True DO_SMIF_STACKING=True DO_SMIF_APBS=False DO_OUTPUT_DX=False DO_OUTPUT_MRC=True DO_OUTPUT_CMAP=False"


############################# CONVERSIONS
# shellcheck disable=SC2086
python3 -W ignore smiffer.py prot $fpdb/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --debug $args_formats

rm -f $fc/1iqj.meta.json
mv $fc/1iqj.stacking.dx $fc/1iqj.stk.dx
mv $fc/1iqj.stacking.mrc $fc/1iqj.stk.mrc
mv $fc/1iqj.stacking.cmap $fc/1iqj.stk.cmap


############################# PACKING
# shellcheck disable=SC2086
python3 -W ignore smiffer.py rna $fpdb/2esj.pdb -o $fp -rxyz 15.708 21.865 -6.397 16.946 --debug $args_smifs

rm -f $fp/2esj.meta.json
mv $fp/2esj.hbacceptors.mrc $fp/2esj.hba.mrc
mv $fp/2esj.hbdonors.mrc    $fp/2esj.hbd.mrc
mv $fp/2esj.hydrophilic.mrc $fp/2esj.phi.mrc
mv $fp/2esj.hydrophobic.mrc $fp/2esj.pho.mrc
mv $fp/2esj.stacking.mrc    $fp/2esj.stk.mrc


############################# UNPACKING
# shellcheck disable=SC2086
python3 -W ignore smiffer.py prot $fpdb/1iqj.pdb -o $fu -rxyz 14.675 4.682 21.475 7.161 --debug $args_smifs

cd $fu
python3 ../../../../smiffer.py pack -i 1iqj.hbacceptors.mrc 1iqj.hbdonors.mrc 1iqj.hydrophilic.mrc 1iqj.hydrophobic.mrc 1iqj.stacking.mrc -o 1iqj.cmap
cd ../../../..

rm -f $fu/1iqj.meta.json $fu/1iqj.*.mrc


############################# FIX CMAP
python3 smiffer.py pack -i $ff/_frames/smiffer_126.hbdonors.cmap $ff/_frames/smiffer_142.hbdonors.cmap $ff/_frames/smiffer_3.hbdonors.cmap $ff/_frames/smiffer_127.hbdonors.cmap $ff/_frames/smiffer_32.hbdonors.cmap $ff/_frames/smiffer_50.hbdonors.cmap -o $ff/hbdonors.issue.cmap
