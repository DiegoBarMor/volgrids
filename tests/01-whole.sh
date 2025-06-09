#!/bin/bash
set -eu

echo
echo ">>> TEST 01: Benchmark, Whole mode"

fapbs="data/apbs"
fpdb="data/pdb"
fout="data/tests/01-whole"
mkdir -p $fout

python3 -W ignore smiffer.py rna  $fpdb/1akx.pdb  -o $fout -a $fapbs/1akx.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/1bg0.pdb  -o $fout -a $fapbs/1bg0.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/1eby.pdb  -o $fout -a $fapbs/1eby.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/1ehe.pdb  -o $fout -a $fapbs/1ehe.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/1h7l.pdb  -o $fout -a $fapbs/1h7l.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/1i9v.pdb  -o $fout -a $fapbs/1i9v.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/1iqj.pdb  -o $fout -a $fapbs/1iqj.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/1ofz.pdb  -o $fout -a $fapbs/1ofz.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/2esj.pdb  -o $fout -a $fapbs/2esj.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/3dd0.pdb  -o $fout -a $fapbs/3dd0.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/3ee4.pdb  -o $fout -a $fapbs/3ee4.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/4f8u.pdb  -o $fout -a $fapbs/4f8u.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/5bjo.pdb  -o $fout -a $fapbs/5bjo.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/5kx9.pdb  -o $fout -a $fapbs/5kx9.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/5m9w.pdb  -o $fout -a $fapbs/5m9w.pqr.dx
python3 -W ignore smiffer.py prot $fpdb/6e9a.pdb  -o $fout -a $fapbs/6e9a.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/6tf3.pdb  -o $fout -a $fapbs/6tf3.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/7oax0.pdb -o $fout -a $fapbs/7oax0.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/7oax1.pdb -o $fout -a $fapbs/7oax1.pqr.dx
python3 -W ignore smiffer.py rna  $fpdb/8eyv.pdb  -o $fout -a $fapbs/8eyv.pqr.dx

rm -f $fout/*.json $fout/*.pdb

names=(1akx 1bg0 1eby 1ehe 1h7l 1i9v 1iqj 1ofz 2esj 3dd0 3ee4 4f8u 5bjo 5kx9 5m9w 6e9a 6tf3 7oax0 7oax1 8eyv)
for name in "${names[@]}"; do
    mapfile -t paths_grids < <(ls "$fout/$name".*) # get the list of files corresponding to the pdb
    python3 smiffer.py pack -i "${paths_grids[@]}" -o "$fout/$name.cmap"  # pack them into a single CMAP file
    rm -f "${paths_grids[@]}" # remove the individual grid files
    cp "$fpdb/$name.pdb" "$fout/$name.pdb"
done
