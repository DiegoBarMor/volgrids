#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 0: Toy systems"

fpdb="testdata/_input/toy_systems"
fout="testdata/smiffer/toy_systems"
tmp_config_ignore_h="$fpdb/ignore_h.config.tmp"
rm -rf $fout; mkdir -p $fout

cp "$fpdb"/*.pdb "$fout/"
cp "$fpdb/ribose_gua.pdb" "$fout/ribose_gua_no_h.pdb"
cp "$fout/peptide.pdb" "$fout/peptide_no_h.pdb" # for visualization purposes

cat > $tmp_config_ignore_h <<- EOM
[SMIFFER]
USE_STRUCTURE_HYDROGENS=False
EOM

python3 run/smiffer.py prot $fout/peptide_no_h.pdb    -o $fout -c $tmp_config_ignore_h
python3 run/smiffer.py prot $fout/peptide.pdb         -o $fout
python3 run/smiffer.py rna  $fout/guanine.pdb         -o $fout
python3 run/smiffer.py rna  $fout/ribose_gua_no_h.pdb -o $fout -c $tmp_config_ignore_h
python3 run/smiffer.py rna  $fout/ribose_gua.pdb      -o $fout
python3 run/smiffer.py rna  $fout/uuu.pdb             -o $fout

python3 run/smiffer.py prot $fout/all_arg.pdb -o $fout
python3 run/smiffer.py prot $fout/all_asn.pdb -o $fout
python3 run/smiffer.py rna  $fout/all_cyt.pdb -o $fout
python3 run/smiffer.py rna  $fout/all_ump.pdb -o $fout

rm -f $tmp_config_ignore_h
