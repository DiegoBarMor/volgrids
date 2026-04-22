#!/bin/bash
set -eu

### folders that should already exist
fframes="testdata/vgtools/_inconsistent_frames"
fpdb_nosolv="testdata/smiffer/pdb_clean"
f_interface="testdata/smiffer/interfaces/prot_rna"
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
fop="$fvgtools/operations"

mkdir -p $fc $fp $fu $ff $fop
rm -f $fc/* $fp/* $fu/* $ff/*.cmap $fop/*.cmap

conf_no_apbs="GRID_FORMAT_OUTPUT=MRC DO_SMIF_STACKING=True DO_SMIF_HBA=True DO_SMIF_HBD=True DO_SMIF_HYDROPHOBIC=True DO_SMIF_HYDROPHILIC=True DO_SMIF_STACKING=True DO_SMIF_APBS=False SAVE_TRIMMING_MASK=False"
conf_just_stacking="DO_SMIF_STACKING=True DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False SAVE_TRIMMING_MASK=False"
conf_dx="$conf_just_stacking GRID_FORMAT_OUTPUT=DX"
conf_mrc="$conf_just_stacking GRID_FORMAT_OUTPUT=MRC"
conf_ccp4="$conf_just_stacking GRID_FORMAT_OUTPUT=CCP4"
conf_cmap="$conf_just_stacking GRID_FORMAT_OUTPUT=CMAP"


############################# CONVERSIONS
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fc -s 4.682 21.475 7.161 14.675 --config "$conf_dx"
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fc -s 4.682 21.475 7.161 14.675 --config "$conf_mrc"
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fc -s 4.682 21.475 7.161 14.675 --config "$conf_ccp4"
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fc -s 4.682 21.475 7.161 14.675 --config "$conf_cmap"

mv $fc/1iqj.stacking.dx   $fc/1iqj.stk.dx
mv $fc/1iqj.stacking.mrc  $fc/1iqj.stk.mrc
mv $fc/1iqj.stacking.ccp4 $fc/1iqj.stk.ccp4
mv $fc/1iqj.stacking.cmap $fc/1iqj.stk.cmap


############################# PACKING
python3 volgrids smiffer rna $fpdb_nosolv/2esj.pdb -o $fp -s 21.865 -6.397 16.946 15.708 --config "$conf_no_apbs"

mv $fp/2esj.hbacceptors.mrc $fp/2esj.hba.mrc
mv $fp/2esj.hbdonors.mrc    $fp/2esj.hbd.mrc
mv $fp/2esj.hydrophilic.mrc $fp/2esj.phi.mrc
mv $fp/2esj.hydrophobic.mrc $fp/2esj.pho.mrc
mv $fp/2esj.stacking.mrc    $fp/2esj.stk.mrc


############################# UNPACKING
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fu -s 4.682 21.475 7.161 14.675 --config "$conf_no_apbs"

cd $fu
python3 ../../../volgrids vgtools pack \
    1iqj.hbacceptors.mrc 1iqj.hbdonors.mrc 1iqj.hydrophilic.mrc 1iqj.hydrophobic.mrc 1iqj.stacking.mrc \
    -o 1iqj.cmap
cd ../../..

rm -f $fu/1iqj.*.mrc


############################# FIX CMAP
python3 volgrids vgtools pack \
    $fframes/smiffer_126.hbdonors.cmap $fframes/smiffer_142.hbdonors.cmap \
    $fframes/smiffer_3.hbdonors.cmap   $fframes/smiffer_127.hbdonors.cmap \
    $fframes/smiffer_32.hbdonors.cmap  $fframes/smiffer_50.hbdonors.cmap  \
    -o $ff/hbdonors.issue.cmap


############################# OPERATIONS
python3 volgrids smiffer prot $f_interface/prot.pdb --config "$conf_just_stacking" -o $fop
python3 volgrids smiffer rna  $f_interface/rna.pdb  --config "$conf_just_stacking" -o $fop
