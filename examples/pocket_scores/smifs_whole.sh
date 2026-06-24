#!/bin/bash
set -eu

### TEMPLATE: change this sh accordingly if using another dataset

fpdb="testdata/smiffer/pdb_clean"
fout="testdata/smiffer/whole"
config="examples/pocket_scores/config.ini"
rm -rf $fout; mkdir -p $fout

python3 volgrids smiffer $fpdb/1akx.pdb  -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/1bg0.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/1eby.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/1ehe.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/1h7l.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/1i9v.pdb  -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/1iqj.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/1ofz.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/2esj.pdb  -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/3dd0.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/3ee4.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/4f8u.pdb  -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/5bjo.pdb  -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/5kx9.pdb  -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/5m9w.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/6e9a.pdb  -o $fout -c $config --pack # prot
python3 volgrids smiffer $fpdb/6tf3.pdb  -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/7oax0.pdb -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/7oax1.pdb -o $fout -c $config --pack # rna
python3 volgrids smiffer $fpdb/8eyv.pdb  -o $fout -c $config --pack # rna
