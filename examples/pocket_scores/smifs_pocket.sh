#!/bin/bash
set -eu

### TEMPLATE: change this sh accordingly if using another dataset

fpdb="testdata/smiffer/pdb-nosolv"
fout="testdata/smiffer/pocket_sphere"
rm -rf $fout; mkdir -p $fout

python3 smiffer.py rna  $fpdb/1akx.pdb  -o $fout -a -s  -2.677  -4.466   -0.020   9.698
python3 smiffer.py prot $fpdb/1bg0.pdb  -o $fout -a -s  21.078  12.898   40.207  11.895
python3 smiffer.py prot $fpdb/1eby.pdb  -o $fout -a -s  13.128  22.854    5.661  10.310
python3 smiffer.py prot $fpdb/1ehe.pdb  -o $fout -a -s  50.266  11.982   20.951  11.555
python3 smiffer.py prot $fpdb/1h7l.pdb  -o $fout -a -s  18.424  16.846   35.302   8.523
python3 smiffer.py rna  $fpdb/1i9v.pdb  -o $fout -a -s  22.853  -8.022   15.465  14.252
python3 smiffer.py prot $fpdb/1iqj.pdb  -o $fout -a -s   4.682  21.475    7.161  14.675
python3 smiffer.py prot $fpdb/1ofz.pdb  -o $fout -a -s  41.092  10.237   74.244   9.563
python3 smiffer.py rna  $fpdb/2esj.pdb  -o $fout -a -s  21.865  -6.397   16.946  15.708
python3 smiffer.py prot $fpdb/3dd0.pdb  -o $fout -a -s  -4.599   3.904   15.646   8.944
python3 smiffer.py prot $fpdb/3ee4.pdb  -o $fout -a -s  28.062  17.653   14.159  12.542
python3 smiffer.py rna  $fpdb/4f8u.pdb  -o $fout -a -s -33.649   0.380   -8.224  12.889
python3 smiffer.py rna  $fpdb/5bjo.pdb  -o $fout -a -s   3.924  10.347  -13.364   6.632
python3 smiffer.py rna  $fpdb/5kx9.pdb  -o $fout -a -s  28.204  38.318   21.295  11.763
python3 smiffer.py prot $fpdb/5m9w.pdb  -o $fout -a -s  58.806  39.475    7.084  11.003
python3 smiffer.py prot $fpdb/6e9a.pdb  -o $fout -a -s  16.802  23.180   17.191  12.875
python3 smiffer.py rna  $fpdb/6tf3.pdb  -o $fout -a -s  21.474  -9.821   18.388  11.750
python3 smiffer.py rna  $fpdb/7oax0.pdb -o $fout -a -s -16.980   7.444   -4.446  12.179
python3 smiffer.py rna  $fpdb/7oax1.pdb -o $fout -a -s  10.243  -5.151   -8.194  14.835
python3 smiffer.py rna  $fpdb/8eyv.pdb  -o $fout -a -s  -1.612  -8.183   18.333  11.998
