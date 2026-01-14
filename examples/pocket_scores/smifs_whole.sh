#!/bin/bash
set -eu

### TEMPLATE: change this sh accordingly if using another dataset

fpdb="testdata/smiffer/pdb-nosolv"
fout="testdata/smiffer/whole"
rm -rf $fout; mkdir -p $fout

python3 smiffer.py rna  $fpdb/1akx.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/1bg0.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/1eby.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/1ehe.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/1h7l.pdb  -o $fout -a
python3 smiffer.py rna  $fpdb/1i9v.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/1iqj.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/1ofz.pdb  -o $fout -a
python3 smiffer.py rna  $fpdb/2esj.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/3dd0.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/3ee4.pdb  -o $fout -a
python3 smiffer.py rna  $fpdb/4f8u.pdb  -o $fout -a
python3 smiffer.py rna  $fpdb/5bjo.pdb  -o $fout -a
python3 smiffer.py rna  $fpdb/5kx9.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/5m9w.pdb  -o $fout -a
python3 smiffer.py prot $fpdb/6e9a.pdb  -o $fout -a
python3 smiffer.py rna  $fpdb/6tf3.pdb  -o $fout -a
python3 smiffer.py rna  $fpdb/7oax0.pdb -o $fout -a
python3 smiffer.py rna  $fpdb/7oax1.pdb -o $fout -a
python3 smiffer.py rna  $fpdb/8eyv.pdb  -o $fout -a
