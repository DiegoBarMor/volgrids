#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 2: Benchmark, Whole mode"

fapbs="testdata/smiffer/apbs"
fpdb_orig="testdata/smiffer/pdb_orig"
fpdb_clean="testdata/smiffer/pdb_clean"
fout="testdata/smiffer/whole"
mkdir -p $fout

names=(1akx 1bg0 6e9a 6tf3 7oax0)
for name in "${names[@]}"; do
    cp "$fpdb_clean/$name.pdb" "$fout/$name.pdb"
done

python3 volgrids smiffer rna  $fpdb_clean/1akx.pdb  -o $fout -a $fapbs/1akx.pdb.mrc
python3 volgrids smiffer prot $fpdb_clean/1bg0.pdb  -o $fout -a $fapbs/1bg0.pdb.mrc
python3 volgrids smiffer prot $fpdb_clean/6e9a.pdb  -o $fout -a $fapbs/6e9a.pdb.mrc
python3 volgrids smiffer rna  $fpdb_clean/7oax0.pdb -o $fout -a $fapbs/7oax0.pdb.mrc


fout="testdata/smiffer/whole-simple_hbonds"
mkdir -p $fout

conf_just_hbond="DO_SMIF_STACKING=false DO_SMIF_HBA=true DO_SMIF_HBD=true DO_SMIF_HYDROPHOBIC=false DO_SMIF_HYDROPHILIC=false DO_SMIF_APBS=false"

names=(1akx 1i9v 2esj 4f8u 5bjo 5kx9 6tf3 7oax0 8eyv)
for name in "${names[@]}"; do
    cp "$fpdb_orig/$name.pdb" "$fout/$name.pdb"
    python3 volgrids smiffer rna  "$fout/$name.pdb" -c "$conf_just_hbond" DO_SIMPLE_HBONDS_RNA=true
    mv "$fout/$name.cmap" "$fout/$name.simple.cmap"
    python3 volgrids smiffer rna  "$fout/$name.pdb" -c "$conf_just_hbond"
done
