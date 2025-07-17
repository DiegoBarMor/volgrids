#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 0: Toy systems"

fpdb="testdata/_input/toy_systems"
fout="testdata/smiffer/toy_systems"
rm -rf $fout; mkdir -p $fout

cp "$fpdb"/*.pdb "$fout/"
cp "$fpdb/ribose_gua.pdb" "$fout/ribose_gua_H.pdb"

python3 smiffer.py prot $fout/peptide.pdb      -o $fout
python3 smiffer.py prot $fout/peptide_H.pdb    -o $fout --protonated
python3 smiffer.py rna  $fout/guanine.pdb      -o $fout
python3 smiffer.py rna  $fout/ribose_gua.pdb   -o $fout
python3 smiffer.py rna  $fout/ribose_gua_H.pdb -o $fout

python3 smiffer.py prot $fout/arg.pdb -o $fout
python3 smiffer.py prot $fout/asn.pdb -o $fout
python3 smiffer.py rna  $fout/cyt.pdb -o $fout
