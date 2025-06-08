#!/bin/bash
set -eu

echo
echo ">>> TEST 00: Benchmark, Pocket Sphere mode"

fapbs="data/apbs"
fpdb="data/pdb"
fout="data/tests/00-pocket_sphere"

mkdir -p $fout
cp $fpdb/*.pdb $fout/

python3 -W ignore smiffer.py rna  $fout/1akx.pdb  -a $fapbs/1akx.pqr.dx  -rxyz  9.698   -2.677  -4.466   -0.020
python3 -W ignore smiffer.py prot $fout/1bg0.pdb  -a $fapbs/1bg0.pqr.dx  -rxyz 11.895   21.078  12.898   40.207
python3 -W ignore smiffer.py prot $fout/1eby.pdb  -a $fapbs/1eby.pqr.dx  -rxyz 10.310   13.128  22.854    5.661
python3 -W ignore smiffer.py prot $fout/1ehe.pdb  -a $fapbs/1ehe.pqr.dx  -rxyz 11.555   50.266  11.982   20.951
python3 -W ignore smiffer.py prot $fout/1h7l.pdb  -a $fapbs/1h7l.pqr.dx  -rxyz  8.523   18.424  16.846   35.302
python3 -W ignore smiffer.py rna  $fout/1i9v.pdb  -a $fapbs/1i9v.pqr.dx  -rxyz 14.252   22.853  -8.022   15.465
python3 -W ignore smiffer.py prot $fout/1iqj.pdb  -a $fapbs/1iqj.pqr.dx  -rxyz 14.675    4.682  21.475    7.161
python3 -W ignore smiffer.py prot $fout/1ofz.pdb  -a $fapbs/1ofz.pqr.dx  -rxyz  9.563   41.092  10.237   74.244
python3 -W ignore smiffer.py rna  $fout/2esj.pdb  -a $fapbs/2esj.pqr.dx  -rxyz 15.708   21.865  -6.397   16.946
python3 -W ignore smiffer.py prot $fout/3dd0.pdb  -a $fapbs/3dd0.pqr.dx  -rxyz  8.944   -4.599   3.904   15.646
python3 -W ignore smiffer.py prot $fout/3ee4.pdb  -a $fapbs/3ee4.pqr.dx  -rxyz 12.542   28.062  17.653   14.159
python3 -W ignore smiffer.py rna  $fout/4f8u.pdb  -a $fapbs/4f8u.pqr.dx  -rxyz 12.889  -33.649   0.380   -8.224
python3 -W ignore smiffer.py rna  $fout/5bjo.pdb  -a $fapbs/5bjo.pqr.dx  -rxyz  6.632    3.924  10.347  -13.364
python3 -W ignore smiffer.py rna  $fout/5kx9.pdb  -a $fapbs/5kx9.pqr.dx  -rxyz 11.763   28.204  38.318   21.295
python3 -W ignore smiffer.py prot $fout/5m9w.pdb  -a $fapbs/5m9w.pqr.dx  -rxyz 11.003   58.806  39.475    7.084
python3 -W ignore smiffer.py prot $fout/6e9a.pdb  -a $fapbs/6e9a.pqr.dx  -rxyz 12.875   16.802  23.180   17.191
python3 -W ignore smiffer.py rna  $fout/6tf3.pdb  -a $fapbs/6tf3.pqr.dx  -rxyz 11.750   21.474  -9.821   18.388
python3 -W ignore smiffer.py rna  $fout/7oax0.pdb -a $fapbs/7oax0.pqr.dx -rxyz 12.179  -16.980   7.444   -4.446
python3 -W ignore smiffer.py rna  $fout/7oax1.pdb -a $fapbs/7oax1.pqr.dx -rxyz 14.835   10.243  -5.151   -8.194
python3 -W ignore smiffer.py rna  $fout/8eyv.pdb  -a $fapbs/8eyv.pqr.dx  -rxyz 11.998   -1.612  -8.183   18.333

rm -f $fout/*.json # remove metadata

names=(1akx 1bg0 1eby 1ehe 1h7l 1i9v 1iqj 1ofz 2esj 3dd0 3ee4 4f8u 5bjo 5kx9 5m9w 6e9a 6tf3 7oax0 7oax1 8eyv)
for name in "${names[@]}"; do
    mapfile -t paths_grids < <(ls "$fout/$name".*) # get the list of files corresponding to the pdb
    python3 smif_tools.py -p "$fout/$name.cmap" "${paths_grids[@]}" # pack them into a single CMAP file
    rm -f "${paths_grids[@]}" # remove the individual grid files
done
