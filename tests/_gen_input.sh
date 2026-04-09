#!/bin/bash
set -eu

### folders that should already exist
fframes="testdata/vgtools/_inconsistent_frames"
fpdb_nosolv="testdata/smiffer/pdb_clean"
fapbs="testdata/smiffer/apbs"

mkdir -p $fapbs

################################################################################
##################################### APBS #####################################
################################################################################

names=(1akx 1bg0 1eby 1ehe 1h7l 1i9v 1iqj 1ofz 2esj 3dd0 3ee4 4f8u 5bjo 5kx9 5m9w 6e9a 6tf3 7oax0 7oax1 8eyv)
for name in "${names[@]}"; do
    path_pdb_orig="$fpdb_nosolv/$name.pdb"
    path_pdb_apbs="$fapbs/$name.pdb"
    cp "$path_pdb_orig" "$path_pdb_apbs"
    python3 volgrids apbs "$path_pdb_apbs" --mrc --verbose
    rm -f "$path_pdb_apbs"
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
GRID_FORMAT_OUTPUT=DX
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM

cat > $tmp_config_mrc <<- EOM
GRID_FORMAT_OUTPUT=MRC
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM

cat > $tmp_config_ccp4 <<- EOM
GRID_FORMAT_OUTPUT=CCP4
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM

cat > $tmp_config_cmap <<- EOM
GRID_FORMAT_OUTPUT=CMAP
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
SAVE_TRIMMING_MASK=False
EOM

cat > $tmp_config_smifs <<- EOM
GRID_FORMAT_OUTPUT=MRC
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
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fc -s 4.682 21.475 7.161 14.675 --config $tmp_config_dx
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fc -s 4.682 21.475 7.161 14.675 --config $tmp_config_mrc
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fc -s 4.682 21.475 7.161 14.675 --config $tmp_config_ccp4
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fc -s 4.682 21.475 7.161 14.675 --config $tmp_config_cmap

mv $fc/1iqj.stacking.dx   $fc/1iqj.stk.dx
mv $fc/1iqj.stacking.mrc  $fc/1iqj.stk.mrc
mv $fc/1iqj.stacking.ccp4 $fc/1iqj.stk.ccp4
mv $fc/1iqj.stacking.cmap $fc/1iqj.stk.cmap


############################# PACKING
# shellcheck disable=SC2086
python3 volgrids smiffer rna $fpdb_nosolv/2esj.pdb -o $fp -s 21.865 -6.397 16.946 15.708 --config $tmp_config_smifs

mv $fp/2esj.hbacceptors.mrc $fp/2esj.hba.mrc
mv $fp/2esj.hbdonors.mrc    $fp/2esj.hbd.mrc
mv $fp/2esj.hydrophilic.mrc $fp/2esj.phi.mrc
mv $fp/2esj.hydrophobic.mrc $fp/2esj.pho.mrc
mv $fp/2esj.stacking.mrc    $fp/2esj.stk.mrc


############################# UNPACKING
# shellcheck disable=SC2086
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fu -s 4.682 21.475 7.161 14.675 --config $tmp_config_smifs

cd $fu
python3 ../../../volgrids vgtools pack \
    -i 1iqj.hbacceptors.mrc 1iqj.hbdonors.mrc 1iqj.hydrophilic.mrc 1iqj.hydrophobic.mrc 1iqj.stacking.mrc \
    -o 1iqj.cmap
cd ../../..

rm -f $fu/1iqj.*.mrc


############################# FIX CMAP
python3 volgrids vgtools pack -i \
    $fframes/smiffer_126.hbdonors.cmap $fframes/smiffer_142.hbdonors.cmap \
    $fframes/smiffer_3.hbdonors.cmap   $fframes/smiffer_127.hbdonors.cmap \
    $fframes/smiffer_32.hbdonors.cmap  $fframes/smiffer_50.hbdonors.cmap  \
    -o $ff/hbdonors.issue.cmap


############################# Cleanup
rm -f $tmp_config_dx $tmp_config_mrc $tmp_config_ccp4 $tmp_config_cmap $tmp_config_smifs
