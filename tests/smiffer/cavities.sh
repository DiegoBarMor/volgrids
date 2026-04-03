#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 5: Cavities"

fpdb="testdata/smiffer/pdb-nosolv"
fout="testdata/smiffer/cavities"
tmp_config="$fout/config.tmp"
rm -rf $fout; mkdir -p $fout

cat > $tmp_config <<- EOM
DO_SMIF_STACKING=False
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False
DO_CAVITIES_FINDER=True
SAVE_TRIMMING_MASK=True
EOM

for name in 1bg0 1eby 1ehe 1h7l 1iqj 1ofz 3dd0 3ee4 5m9w 6e9a ; do
    python3 volgrids smiffer prot  $fpdb/$name.pdb -o $fout --config $tmp_config
    cp $fpdb/$name.pdb $fout/
done
for name in 1akx 1i9v 2esj 4f8u 5bjo 5kx9 6tf3 7oax0 7oax1 8eyv; do
    python3 volgrids smiffer rna  $fpdb/$name.pdb -o $fout --config $tmp_config
    cp $fpdb/$name.pdb $fout/
done

rm -f $tmp_config
