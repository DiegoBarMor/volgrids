#!/bin/bash
set -eu

fdata="testdata/_input"
fpdb_orig="$fdata/pdb-original"
fpdb_nosolv="$fdata/pdb-nosolv"
fapbs="$fdata/apbs"
fframes="$fdata/_example_frames"



################################################################################
##################################### APBS #####################################
################################################################################
rm -rf $fpdb_nosolv $fapbs
mkdir -p $fpdb_nosolv $fapbs

tmp_pyscript=$fdata/remove_solvent.tmp.py
cat > $tmp_pyscript <<- EOM
import sys
import MDAnalysis as mda
from pathlib import Path
folder_in = Path(sys.argv[1]); folder_out = Path(sys.argv[2])
benchmark_prots = ["1bg0","1eby","1ehe","1h7l","1iqj","1ofz","3dd0","3ee4", "5m9w", "6e9a"]
benchmark_rnas  = ["1akx","1i9v","2esj","4f8u","5bjo","5kx9","6tf3","7oax0","7oax1","8eyv"]
def remove_solv(path_in, path_out, selection): mda.Universe(path_in).select_atoms(selection).write(path_out)
for name in benchmark_prots: remove_solv(folder_in / f"{name}.pdb", folder_out / f"{name}.pdb", "protein")
for name in benchmark_rnas:  remove_solv(folder_in / f"{name}.pdb", folder_out / f"{name}.pdb", "nucleic")
EOM
echo "Removing solvent from PDB files..."
python3 -W ignore $tmp_pyscript $fpdb_orig $fpdb_nosolv
rm -f $tmp_pyscript

names=(1akx 1bg0 1eby 1ehe 1h7l 1i9v 1iqj 1ofz 2esj 3dd0 3ee4 4f8u 5bjo 5kx9 5m9w 6e9a 6tf3 7oax0 7oax1 8eyv)
for name in "${names[@]}"; do
    bash run_apbs.sh "$fpdb_nosolv/$name.pdb" $fapbs
done



################################################################################
#################################### VGTOOLS ###################################
################################################################################

fvgtools="testdata/vgtools"
fc="$fvgtools/converting"
fp="$fvgtools/packing"
fu="$fvgtools/unpacking"
ff="$fvgtools/fix_cmap"

mkdir -p $fc $fp $fu $ff
rm -f $fc/* $fp/* $fu/*  $ff/*.cmap

tmp_config_dx=$fvgtools/dx.config
tmp_config_mrc=$fvgtools/mrc.config
tmp_config_ccp4=$fvgtools/ccp4.config
tmp_config_cmap=$fvgtools/cmap.config
tmp_config_smifs=$fvgtools/smifs.config

cat > $tmp_config_dx <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=GridFormat.DX
[SMIFFER]
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM

cat > $tmp_config_mrc <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=GridFormat.MRC
[SMIFFER]
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM

cat > $tmp_config_ccp4 <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=GridFormat.CCP4
[SMIFFER]
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM

cat > $tmp_config_cmap <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=GridFormat.CMAP
[SMIFFER]
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM

cat > $tmp_config_smifs <<- EOM
[VOLGRIDS]
OUTPUT_FORMAT=GridFormat.MRC
[SMIFFER]
DO_SMIF_HBA=True
DO_SMIF_HBD=True
DO_SMIF_HYDROPHOBIC=True
DO_SMIF_HYDROPHILIC=True
DO_SMIF_STACKING=True
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM


############################# CONVERSIONS
# shellcheck disable=SC2086
python3 -W ignore smiffer.py prot $fpdb_nosolv/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_dx
python3 -W ignore smiffer.py prot $fpdb_nosolv/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_mrc
python3 -W ignore smiffer.py prot $fpdb_nosolv/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_ccp4
python3 -W ignore smiffer.py prot $fpdb_nosolv/1iqj.pdb -o $fc -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_cmap

mv $fc/1iqj.stacking.dx   $fc/1iqj.stk.dx
mv $fc/1iqj.stacking.mrc  $fc/1iqj.stk.mrc
mv $fc/1iqj.stacking.ccp4 $fc/1iqj.stk.ccp4
mv $fc/1iqj.stacking.cmap $fc/1iqj.stk.cmap


############################# PACKING
# shellcheck disable=SC2086
python3 -W ignore smiffer.py rna $fpdb_nosolv/2esj.pdb -o $fp -rxyz 15.708 21.865 -6.397 16.946 --config $tmp_config_smifs

mv $fp/2esj.hbacceptors.mrc $fp/2esj.hba.mrc
mv $fp/2esj.hbdonors.mrc    $fp/2esj.hbd.mrc
mv $fp/2esj.hydrophilic.mrc $fp/2esj.phi.mrc
mv $fp/2esj.hydrophobic.mrc $fp/2esj.pho.mrc
mv $fp/2esj.stacking.mrc    $fp/2esj.stk.mrc


############################# UNPACKING
# shellcheck disable=SC2086
python3 -W ignore smiffer.py prot $fpdb_nosolv/1iqj.pdb -o $fu -rxyz 14.675 4.682 21.475 7.161 --config $tmp_config_smifs

cd $fu
python3 ../../../vgtools.py pack \
    -i 1iqj.hbacceptors.mrc 1iqj.hbdonors.mrc 1iqj.hydrophilic.mrc 1iqj.hydrophobic.mrc 1iqj.stacking.mrc \
    -o 1iqj.cmap
cd ../../..

rm -f $fu/1iqj.*.mrc


############################# FIX CMAP
python3 vgtools.py pack -i \
    $fframes/smiffer_126.hbdonors.cmap $fframes/smiffer_142.hbdonors.cmap \
    $fframes/smiffer_3.hbdonors.cmap   $fframes/smiffer_127.hbdonors.cmap \
    $fframes/smiffer_32.hbdonors.cmap  $fframes/smiffer_50.hbdonors.cmap  \
    -o $ff/hbdonors.issue.cmap


############################# Cleanup
rm -f $tmp_config_dx $tmp_config_mrc $tmp_config_ccp4 $tmp_config_cmap $tmp_config_smifs
