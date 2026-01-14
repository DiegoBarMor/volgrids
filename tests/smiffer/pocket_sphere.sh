#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 1: Benchmark, Pocket Sphere mode"

fapbs="testdata/smiffer/apbs"
fpdb="testdata/smiffer/pdb-nosolv"
fout="testdata/smiffer/pocket_sphere"
rm -rf $fout; mkdir -p $fout

names=(1akx 1bg0 1eby 1ehe 1h7l 1i9v 1iqj 1ofz 2esj 3dd0 3ee4 4f8u 5bjo 5kx9 5m9w 6e9a 6tf3 7oax0 7oax1 8eyv)
for name in "${names[@]}"; do
    cp "$fpdb/$name.pdb" "$fout/$name.pdb"
done

python3 volgrids smiffer rna  $fpdb/1akx.pdb  -o $fout -a $fapbs/1akx.pdb.mrc  -s  -2.677  -4.466   -0.020   9.698
python3 volgrids smiffer prot $fpdb/1bg0.pdb  -o $fout -a $fapbs/1bg0.pdb.mrc  -s  21.078  12.898   40.207  11.895
python3 volgrids smiffer prot $fpdb/1eby.pdb  -o $fout -a $fapbs/1eby.pdb.mrc  -s  13.128  22.854    5.661  10.310
python3 volgrids smiffer prot $fpdb/1ehe.pdb  -o $fout -a $fapbs/1ehe.pdb.mrc  -s  50.266  11.982   20.951  11.555
python3 volgrids smiffer prot $fpdb/1h7l.pdb  -o $fout -a $fapbs/1h7l.pdb.mrc  -s  18.424  16.846   35.302   8.523
python3 volgrids smiffer rna  $fpdb/1i9v.pdb  -o $fout -a $fapbs/1i9v.pdb.mrc  -s  22.853  -8.022   15.465  14.252
python3 volgrids smiffer prot $fpdb/1iqj.pdb  -o $fout -a $fapbs/1iqj.pdb.mrc  -s   4.682  21.475    7.161  14.675
python3 volgrids smiffer prot $fpdb/1ofz.pdb  -o $fout -a $fapbs/1ofz.pdb.mrc  -s  41.092  10.237   74.244   9.563
python3 volgrids smiffer rna  $fpdb/2esj.pdb  -o $fout -a $fapbs/2esj.pdb.mrc  -s  21.865  -6.397   16.946  15.708
python3 volgrids smiffer prot $fpdb/3dd0.pdb  -o $fout -a $fapbs/3dd0.pdb.mrc  -s  -4.599   3.904   15.646   8.944
python3 volgrids smiffer prot $fpdb/3ee4.pdb  -o $fout -a $fapbs/3ee4.pdb.mrc  -s  28.062  17.653   14.159  12.542
python3 volgrids smiffer rna  $fpdb/4f8u.pdb  -o $fout -a $fapbs/4f8u.pdb.mrc  -s -33.649   0.380   -8.224  12.889
python3 volgrids smiffer rna  $fpdb/5bjo.pdb  -o $fout -a $fapbs/5bjo.pdb.mrc  -s   3.924  10.347  -13.364   6.632
python3 volgrids smiffer rna  $fpdb/5kx9.pdb  -o $fout -a $fapbs/5kx9.pdb.mrc  -s  28.204  38.318   21.295  11.763
python3 volgrids smiffer prot $fpdb/5m9w.pdb  -o $fout -a $fapbs/5m9w.pdb.mrc  -s  58.806  39.475    7.084  11.003
python3 volgrids smiffer prot $fpdb/6e9a.pdb  -o $fout -a $fapbs/6e9a.pdb.mrc  -s  16.802  23.180   17.191  12.875
python3 volgrids smiffer rna  $fpdb/6tf3.pdb  -o $fout -a $fapbs/6tf3.pdb.mrc  -s  21.474  -9.821   18.388  11.750
python3 volgrids smiffer rna  $fpdb/7oax0.pdb -o $fout -a $fapbs/7oax0.pdb.mrc -s -16.980   7.444   -4.446  12.179
python3 volgrids smiffer rna  $fpdb/7oax1.pdb -o $fout -a $fapbs/7oax1.pdb.mrc -s  10.243  -5.151   -8.194  14.835
python3 volgrids smiffer rna  $fpdb/8eyv.pdb  -o $fout -a $fapbs/8eyv.pdb.mrc  -s  -1.612  -8.183   18.333  11.998

# cp "$fpdb/2g5k.pdb" "$fout/"
# cp "$fpdb/2o3v.pdb" "$fout/"
# python3 volgrids smiffer rna $fpdb/2g5k.pdb -o $fout -s 3.476 27.8915 66.1055 15.0 -a $fapbs/2g5k.pdb.mrc
# python3 volgrids smiffer rna $fpdb/2o3v.pdb -o $fout -s 0.478 -7.3925 45.962  15.0 -a $fapbs/2o3v.pdb.mrc
