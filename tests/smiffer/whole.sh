#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 2: Benchmark, Whole mode"

fapbs="testdata/smiffer/apbs"
fpdb="testdata/smiffer/pdb-nosolv"
fout="testdata/smiffer/whole"
rm -rf $fout; mkdir -p $fout

names=(1akx 1bg0 1eby 1ehe 1h7l 1i9v 1iqj 1ofz 2esj 3dd0 3ee4 4f8u 5bjo 5kx9 5m9w 6e9a 6tf3 7oax0 7oax1 8eyv)
for name in "${names[@]}"; do
    cp "$fpdb/$name.pdb" "$fout/$name.pdb"
done

python3 run/smiffer.py rna  $fpdb/1akx.pdb  -o $fout -a $fapbs/1akx.pdb.mrc
python3 run/smiffer.py prot $fpdb/1bg0.pdb  -o $fout -a $fapbs/1bg0.pdb.mrc
python3 run/smiffer.py prot $fpdb/6e9a.pdb  -o $fout -a $fapbs/6e9a.pdb.mrc
python3 run/smiffer.py rna  $fpdb/7oax0.pdb -o $fout -a $fapbs/7oax0.pdb.mrc
