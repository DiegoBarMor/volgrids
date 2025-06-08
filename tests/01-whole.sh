#!/bin/bash
set -eu

echo
echo ">>> TEST 01: Benchmark, Whole mode"

fapbs="data/apbs"
fpdb="data/pdb"
fout="data/tests/01-whole"

mkdir -p $fout
cp $fpdb/*.pdb $fout/

python3 -W ignore smiffer.py rna  $fout/1akx.pdb  -a $fapbs/1akx.pqr.dx
python3 -W ignore smiffer.py prot $fout/1bg0.pdb  -a $fapbs/1bg0.pqr.dx
python3 -W ignore smiffer.py prot $fout/1eby.pdb  -a $fapbs/1eby.pqr.dx
python3 -W ignore smiffer.py prot $fout/1ehe.pdb  -a $fapbs/1ehe.pqr.dx
python3 -W ignore smiffer.py prot $fout/1h7l.pdb  -a $fapbs/1h7l.pqr.dx
python3 -W ignore smiffer.py rna  $fout/1i9v.pdb  -a $fapbs/1i9v.pqr.dx
python3 -W ignore smiffer.py prot $fout/1iqj.pdb  -a $fapbs/1iqj.pqr.dx
python3 -W ignore smiffer.py prot $fout/1ofz.pdb  -a $fapbs/1ofz.pqr.dx
python3 -W ignore smiffer.py rna  $fout/2esj.pdb  -a $fapbs/2esj.pqr.dx
python3 -W ignore smiffer.py prot $fout/3dd0.pdb  -a $fapbs/3dd0.pqr.dx
python3 -W ignore smiffer.py prot $fout/3ee4.pdb  -a $fapbs/3ee4.pqr.dx
python3 -W ignore smiffer.py rna  $fout/4f8u.pdb  -a $fapbs/4f8u.pqr.dx
python3 -W ignore smiffer.py rna  $fout/5bjo.pdb  -a $fapbs/5bjo.pqr.dx
python3 -W ignore smiffer.py rna  $fout/5kx9.pdb  -a $fapbs/5kx9.pqr.dx
python3 -W ignore smiffer.py prot $fout/5m9w.pdb  -a $fapbs/5m9w.pqr.dx
python3 -W ignore smiffer.py prot $fout/6e9a.pdb  -a $fapbs/6e9a.pqr.dx
python3 -W ignore smiffer.py rna  $fout/6tf3.pdb  -a $fapbs/6tf3.pqr.dx
python3 -W ignore smiffer.py rna  $fout/7oax0.pdb -a $fapbs/7oax0.pqr.dx
python3 -W ignore smiffer.py rna  $fout/7oax1.pdb -a $fapbs/7oax1.pqr.dx
python3 -W ignore smiffer.py rna  $fout/8eyv.pdb  -a $fapbs/8eyv.pqr.dx

rm -f $fout/*.json # remove metadata

names=(1akx 1bg0 1eby 1ehe 1h7l 1i9v 1iqj 1ofz 2esj 3dd0 3ee4 4f8u 5bjo 5kx9 5m9w 6e9a 6tf3 7oax0 7oax1 8eyv)
for name in "${names[@]}"; do
    mapfile -t paths_grids < <(ls "$fout/$name".*) # get the list of files corresponding to the pdb
    python3 smif_tools.py -p "$fout/$name.cmap" "${paths_grids[@]}" # pack them into a single CMAP file
    rm -f "${paths_grids[@]}" # remove the individual grid files
done
