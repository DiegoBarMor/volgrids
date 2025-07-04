#!/bin/bash
set -eu

fpdb="testdata/_input/pdb"
ftests="testdata/vgtools"
fc="$ftests/converting"
fp="$ftests/packing"
fu="$ftests/unpacking"
ff="$ftests/fix_cmap"

mkdir -p $fc $fp $fu $ff
rm -f $fc/* $fp/* $fu/*  $ff/*.cmap

tmp_config_dx=$ftests/dx.config
tmp_config_mrc=$ftests/mrc.config
tmp_config_ccp4=$ftests/ccp4.config
tmp_config_cmap=$ftests/cmap.config
tmp_config_smifs=$ftests/smifs.config

cat > $tmp_config_dx <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=DX
[SMIFFER]
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
EOM

cat > $tmp_config_mrc <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=MRC
[SMIFFER]
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
EOM

cat > $tmp_config_ccp4 <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=CCP4
[SMIFFER]
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
EOM

cat > $tmp_config_cmap <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=CMAP
[SMIFFER]
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
EOM

cat > $tmp_config_smifs <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=MRC
[SMIFFER]
DO_SMIF_HBA=True
DO_SMIF_HBD=True
DO_SMIF_HYDROPHOBIC=True
DO_SMIF_HYDROPHILIC=True
DO_SMIF_STACKING=True
DO_SMIF_APBS=False
EOM


############################# CONVERSIONS
# shellcheck disable=SC2086
python3 -W ignore smiffer.py prot $fpdb/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_dx
python3 -W ignore smiffer.py prot $fpdb/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_mrc
python3 -W ignore smiffer.py prot $fpdb/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_ccp4
python3 -W ignore smiffer.py prot $fpdb/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_cmap

mv $fc/1iqj.stacking.dx $fc/1iqj.stk.dx
mv $fc/1iqj.stacking.mrc $fc/1iqj.stk.mrc
mv $fc/1iqj.stacking.ccp4 $fc/1iqj.stk.ccp4
mv $fc/1iqj.stacking.cmap $fc/1iqj.stk.cmap


############################# PACKING
# shellcheck disable=SC2086
python3 -W ignore smiffer.py rna $fpdb/2esj.pdb -o $fp -rxyz 15.708 21.865 -6.397 16.946 --config $tmp_config_smifs

mv $fp/2esj.hbacceptors.mrc $fp/2esj.hba.mrc
mv $fp/2esj.hbdonors.mrc    $fp/2esj.hbd.mrc
mv $fp/2esj.hydrophilic.mrc $fp/2esj.phi.mrc
mv $fp/2esj.hydrophobic.mrc $fp/2esj.pho.mrc
mv $fp/2esj.stacking.mrc    $fp/2esj.stk.mrc


############################# UNPACKING
# shellcheck disable=SC2086
python3 -W ignore smiffer.py prot $fpdb/1iqj.pdb -o $fu -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_smifs

cd $fu
python3 ../../../vgtools.py pack -i 1iqj.hbacceptors.mrc 1iqj.hbdonors.mrc 1iqj.hydrophilic.mrc 1iqj.hydrophobic.mrc 1iqj.stacking.mrc -o 1iqj.cmap
cd ../../..

rm -f $fu/1iqj.*.mrc


############################# FIX CMAP
python3 vgtools.py pack -i $ff/_frames/smiffer_126.hbdonors.cmap $ff/_frames/smiffer_142.hbdonors.cmap $ff/_frames/smiffer_3.hbdonors.cmap $ff/_frames/smiffer_127.hbdonors.cmap $ff/_frames/smiffer_32.hbdonors.cmap $ff/_frames/smiffer_50.hbdonors.cmap -o $ff/hbdonors.issue.cmap


############################# Cleanup
rm -f $tmp_config_dx $tmp_config_mrc $tmp_config_ccp4 $tmp_config_cmap $tmp_config_smifs
