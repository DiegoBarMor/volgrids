#!/bin/bash
set -eu

fapbs="data/apbs"
fpdb="data/pdb"
fout="data/tests/00-pocket_sphere"

python3 -W ignore smiffer.py -i $fpdb/1akx.pdb  -a $fapbs/1akx.pqr.dx  -o $fout -n -r 9.698  -x -2.677  -y -4.466 -z -0.02
python3 -W ignore smiffer.py -i $fpdb/1bg0.pdb  -a $fapbs/1bg0.pqr.dx  -o $fout    -r 11.895 -x 21.078  -y 12.898 -z 40.207
python3 -W ignore smiffer.py -i $fpdb/1eby.pdb  -a $fapbs/1eby.pqr.dx  -o $fout    -r 10.31  -x 13.128  -y 22.854 -z 5.661
python3 -W ignore smiffer.py -i $fpdb/1ehe.pdb  -a $fapbs/1ehe.pqr.dx  -o $fout    -r 11.555 -x 50.266  -y 11.982 -z 20.951
python3 -W ignore smiffer.py -i $fpdb/1h7l.pdb  -a $fapbs/1h7l.pqr.dx  -o $fout    -r 8.523  -x 18.424  -y 16.846 -z 35.302
python3 -W ignore smiffer.py -i $fpdb/1i9v.pdb  -a $fapbs/1i9v.pqr.dx  -o $fout -n -r 14.252 -x 22.853  -y -8.022 -z 15.465
python3 -W ignore smiffer.py -i $fpdb/1iqj.pdb  -a $fapbs/1iqj.pqr.dx  -o $fout    -r 14.675 -x 4.682   -y 21.475 -z 7.161
python3 -W ignore smiffer.py -i $fpdb/1ofz.pdb  -a $fapbs/1ofz.pqr.dx  -o $fout    -r 9.563  -x 41.092  -y 10.237 -z 74.244
python3 -W ignore smiffer.py -i $fpdb/2esj.pdb  -a $fapbs/2esj.pqr.dx  -o $fout -n -r 15.708 -x 21.865  -y -6.397 -z 16.946
python3 -W ignore smiffer.py -i $fpdb/3dd0.pdb  -a $fapbs/3dd0.pqr.dx  -o $fout    -r 8.944  -x -4.599  -y 3.904  -z 15.646
python3 -W ignore smiffer.py -i $fpdb/3ee4.pdb  -a $fapbs/3ee4.pqr.dx  -o $fout    -r 12.542 -x 28.062  -y 17.653 -z 14.159
python3 -W ignore smiffer.py -i $fpdb/4f8u.pdb  -a $fapbs/4f8u.pqr.dx  -o $fout -n -r 12.889 -x -33.649 -y 0.38   -z -8.224
python3 -W ignore smiffer.py -i $fpdb/5bjo.pdb  -a $fapbs/5bjo.pqr.dx  -o $fout -n -r 6.632  -x 3.924   -y 10.347 -z -13.364
python3 -W ignore smiffer.py -i $fpdb/5kx9.pdb  -a $fapbs/5kx9.pqr.dx  -o $fout -n -r 11.763 -x 28.204  -y 38.318 -z 21.295
python3 -W ignore smiffer.py -i $fpdb/5m9w.pdb  -a $fapbs/5m9w.pqr.dx  -o $fout    -r 11.003 -x 58.806  -y 39.475 -z 7.084
python3 -W ignore smiffer.py -i $fpdb/6e9a.pdb  -a $fapbs/6e9a.pqr.dx  -o $fout    -r 12.875 -x 16.802  -y 23.18  -z 17.191
python3 -W ignore smiffer.py -i $fpdb/6tf3.pdb  -a $fapbs/6tf3.pqr.dx  -o $fout -n -r 11.75  -x 21.474  -y -9.821 -z 18.388
python3 -W ignore smiffer.py -i $fpdb/7oax0.pdb -a $fapbs/7oax0.pqr.dx -o $fout -n -r 12.179 -x -16.98  -y 7.444  -z -4.446
python3 -W ignore smiffer.py -i $fpdb/7oax1.pdb -a $fapbs/7oax1.pqr.dx -o $fout -n -r 14.835 -x 10.243  -y -5.151 -z -8.194
python3 -W ignore smiffer.py -i $fpdb/8eyv.pdb  -a $fapbs/8eyv.pqr.dx  -o $fout -n -r 11.998 -x -1.612  -y -8.183 -z 18.333
