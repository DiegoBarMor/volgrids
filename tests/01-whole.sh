#!/bin/bash
set -eu

echo
echo ">>> TEST 01: Benchmark, Whole mode"

fapbs="data/apbs"
fpdb="data/pdb"
fout="data/tests/01-whole"

python3 -W ignore smiffer.py -i $fpdb/1akx.pdb  -a $fapbs/1akx.pqr.dx  -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/1bg0.pdb  -a $fapbs/1bg0.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/1eby.pdb  -a $fapbs/1eby.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/1ehe.pdb  -a $fapbs/1ehe.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/1h7l.pdb  -a $fapbs/1h7l.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/1i9v.pdb  -a $fapbs/1i9v.pqr.dx  -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/1iqj.pdb  -a $fapbs/1iqj.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/1ofz.pdb  -a $fapbs/1ofz.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/2esj.pdb  -a $fapbs/2esj.pqr.dx  -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/3dd0.pdb  -a $fapbs/3dd0.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/3ee4.pdb  -a $fapbs/3ee4.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/4f8u.pdb  -a $fapbs/4f8u.pqr.dx  -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/5bjo.pdb  -a $fapbs/5bjo.pqr.dx  -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/5kx9.pdb  -a $fapbs/5kx9.pqr.dx  -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/5m9w.pdb  -a $fapbs/5m9w.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/6e9a.pdb  -a $fapbs/6e9a.pqr.dx  -o $fout    -w
python3 -W ignore smiffer.py -i $fpdb/6tf3.pdb  -a $fapbs/6tf3.pqr.dx  -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/7oax0.pdb -a $fapbs/7oax0.pqr.dx -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/7oax1.pdb -a $fapbs/7oax1.pqr.dx -o $fout -n -w
python3 -W ignore smiffer.py -i $fpdb/8eyv.pdb  -a $fapbs/8eyv.pqr.dx  -o $fout -n -w
