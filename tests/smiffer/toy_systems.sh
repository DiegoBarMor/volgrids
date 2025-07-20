#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 0: Toy systems"

fpdb="testdata/_input/toy_systems"
fout="testdata/smiffer/toy_systems"
rm -rf $fout; mkdir -p $fout

cp "$fpdb"/*.pdb "$fout/"
cp "$fpdb/ribose_gua.pdb" "$fout/ribose_gua_H.pdb"
cp "$fout/peptide_H.pdb" "$fout/peptide.pdb"

python3 run/smiffer.py prot $fout/peptide.pdb      -o $fout
python3 run/smiffer.py prot $fout/peptide_H.pdb    -o $fout --protonated
python3 run/smiffer.py rna  $fout/guanine.pdb      -o $fout
python3 run/smiffer.py rna  $fout/ribose_gua.pdb   -o $fout
python3 run/smiffer.py rna  $fout/ribose_gua_H.pdb -o $fout --protonated
python3 run/smiffer.py rna  $fout/uuu.pdb          -o $fout

python3 run/smiffer.py prot $fout/all_arg.pdb -o $fout
python3 run/smiffer.py prot $fout/all_asn.pdb -o $fout
python3 run/smiffer.py rna  $fout/all_cyt.pdb -o $fout
python3 run/smiffer.py rna  $fout/all_ump.pdb -o $fout
