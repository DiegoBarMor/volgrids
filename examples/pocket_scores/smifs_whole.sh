#!/bin/bash
set -eu

### TEMPLATE: change this sh accordingly if using another dataset

fpdb="testdata/smiffer/pdb-nosolv"
fout="testdata/smiffer/whole"
config="examples/pocket_scores/config.ini"
rm -rf $fout; mkdir -p $fout

python3 volgrids smiffer rna  $fpdb/1akx.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/1bg0.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/1eby.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/1ehe.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/1h7l.pdb  -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/1i9v.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/1iqj.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/1ofz.pdb  -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/2esj.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/3dd0.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/3ee4.pdb  -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/4f8u.pdb  -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/5bjo.pdb  -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/5kx9.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/5m9w.pdb  -o $fout -c $config
python3 volgrids smiffer prot $fpdb/6e9a.pdb  -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/6tf3.pdb  -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/7oax0.pdb -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/7oax1.pdb -o $fout -c $config
python3 volgrids smiffer rna  $fpdb/8eyv.pdb  -o $fout -c $config
