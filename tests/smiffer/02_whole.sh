#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 2: Benchmark, Whole mode"

fapbs="testdata/smiffer/apbs"
fpdb_clean="testdata/smiffer/pdb_clean"
fout="testdata/smiffer/whole"
mkdir -p $fout

names=(1akx 1bg0 6e9a 6tf3 7oax0)
for name in "${names[@]}"; do
    cp "$fpdb_clean/$name.pdb" "$fout/$name.pdb"
done

python3 volgrids smiffer rna  $fpdb_clean/1akx.pdb  -o $fout -a $fapbs/1akx.pdb.mrc
python3 volgrids smiffer prot $fpdb_clean/1bg0.pdb  -o $fout -a $fapbs/1bg0.pdb.mrc
python3 volgrids smiffer rna  $fpdb_clean/2esj.pdb  -o $fout -a $fapbs/2esj.pdb.mrc
python3 volgrids smiffer prot $fpdb_clean/6e9a.pdb  -o $fout -a $fapbs/6e9a.pdb.mrc
python3 volgrids smiffer rna  $fpdb_clean/7oax0.pdb -o $fout -a $fapbs/7oax0.pdb.mrc
