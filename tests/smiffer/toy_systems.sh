#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 0: Toy systems"

fout="testdata/smiffer/toy_systems"
tmp_config_ignore_h="$fout/ignore_h.config.tmp"
tmp_config_no_apbs="$fout/no_apbs.config.tmp"
tmp_config_equilateral="$fout/equilateral.config.tmp"

cat > $tmp_config_no_apbs <<- EOM
DO_SMIF_APBS=False
EOM

cat > $tmp_config_ignore_h <<- EOM
DO_SMIF_APBS=False
USE_STRUCTURE_HYDROGENS=False
EOM

cat > $tmp_config_equilateral <<- EOM
DO_SMIF_APBS=False
ENSURE_EQUILATERAL=True
EOM

python3 volgrids smiffer prot $fout/peptide_no_h.pdb    -o $fout -c $tmp_config_ignore_h
python3 volgrids smiffer prot $fout/peptide.pdb         -o $fout -c $tmp_config_no_apbs
python3 volgrids smiffer rna  $fout/guanine.pdb         -o $fout -c $tmp_config_no_apbs
python3 volgrids smiffer rna  $fout/ribose_gua_no_h.pdb -o $fout -c $tmp_config_ignore_h
python3 volgrids smiffer rna  $fout/ribose_gua.pdb      -o $fout -c $tmp_config_no_apbs
python3 volgrids smiffer rna  $fout/uuu.pdb             -o $fout -c $tmp_config_no_apbs

python3 volgrids smiffer prot $fout/all_arg.pdb -o $fout -c $tmp_config_no_apbs
python3 volgrids smiffer prot $fout/all_asn.pdb -o $fout -c $tmp_config_no_apbs
python3 volgrids smiffer rna  $fout/all_cyt.pdb -o $fout -c $tmp_config_no_apbs
python3 volgrids smiffer rna  $fout/all_ump.pdb -o $fout -c $tmp_config_no_apbs

cp $fout/all_cyt.pdb $fout/all_cyt_equilateral.pdb
python3 volgrids smiffer rna $fout/all_cyt_equilateral.pdb -o $fout -c $tmp_config_equilateral
rm -f $fout/all_cyt_equilateral.pdb

rm -f $tmp_config_ignore_h $tmp_config_no_apbs $tmp_config_equilateral
