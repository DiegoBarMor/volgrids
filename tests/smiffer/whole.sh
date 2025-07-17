#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 2: Benchmark, Whole mode"

fapbs="testdata/_input/apbs"
fpdb="testdata/_input/pdb-nosolv"
fout="testdata/smiffer/whole"
rm -rf $fout; mkdir -p $fout

names=(1akx 1bg0 1eby 1ehe 1h7l 1i9v 1iqj 1ofz 2esj 3dd0 3ee4 4f8u 5bjo 5kx9 5m9w 6e9a 6tf3 7oax0 7oax1 8eyv)
for name in "${names[@]}"; do
    cp "$fpdb/$name.pdb" "$fout/$name.pdb"
done

python3 run/smiffer.py rna  $fpdb/1akx.pdb  -o $fout -a $fapbs/1akx.pdb.mrc
python3 run/smiffer.py prot $fpdb/1bg0.pdb  -o $fout -a $fapbs/1bg0.pdb.mrc
python3 run/smiffer.py prot $fpdb/1eby.pdb  -o $fout -a $fapbs/1eby.pdb.mrc
python3 run/smiffer.py prot $fpdb/1ehe.pdb  -o $fout -a $fapbs/1ehe.pdb.mrc
python3 run/smiffer.py prot $fpdb/1h7l.pdb  -o $fout -a $fapbs/1h7l.pdb.mrc
python3 run/smiffer.py rna  $fpdb/1i9v.pdb  -o $fout -a $fapbs/1i9v.pdb.mrc
python3 run/smiffer.py prot $fpdb/1iqj.pdb  -o $fout -a $fapbs/1iqj.pdb.mrc
python3 run/smiffer.py prot $fpdb/1ofz.pdb  -o $fout -a $fapbs/1ofz.pdb.mrc
python3 run/smiffer.py rna  $fpdb/2esj.pdb  -o $fout -a $fapbs/2esj.pdb.mrc
python3 run/smiffer.py prot $fpdb/3dd0.pdb  -o $fout -a $fapbs/3dd0.pdb.mrc
python3 run/smiffer.py prot $fpdb/3ee4.pdb  -o $fout -a $fapbs/3ee4.pdb.mrc
python3 run/smiffer.py rna  $fpdb/4f8u.pdb  -o $fout -a $fapbs/4f8u.pdb.mrc
python3 run/smiffer.py rna  $fpdb/5bjo.pdb  -o $fout -a $fapbs/5bjo.pdb.mrc
python3 run/smiffer.py rna  $fpdb/5kx9.pdb  -o $fout -a $fapbs/5kx9.pdb.mrc
python3 run/smiffer.py prot $fpdb/5m9w.pdb  -o $fout -a $fapbs/5m9w.pdb.mrc
python3 run/smiffer.py prot $fpdb/6e9a.pdb  -o $fout -a $fapbs/6e9a.pdb.mrc
python3 run/smiffer.py rna  $fpdb/6tf3.pdb  -o $fout -a $fapbs/6tf3.pdb.mrc
python3 run/smiffer.py rna  $fpdb/7oax0.pdb -o $fout -a $fapbs/7oax0.pdb.mrc
python3 run/smiffer.py rna  $fpdb/7oax1.pdb -o $fout -a $fapbs/7oax1.pdb.mrc
python3 run/smiffer.py rna  $fpdb/8eyv.pdb  -o $fout -a $fapbs/8eyv.pdb.mrc
