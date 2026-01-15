#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 0: Toy systems"

fpdb="testdata/_raw_input/toy_systems"
fout="testdata/smiffer/toy_systems"
tmp_config_ignore_h="$fpdb/ignore_h.config.tmp"
tmp_config_no_apbs="$fpdb/no_apbs.config.tmp"
rm -rf $fout; mkdir -p $fout

cp "$fpdb"/*.pdb "$fout/"
cp "$fpdb/ribose_gua.pdb" "$fout/ribose_gua_no_h.pdb"
cp "$fout/peptide.pdb" "$fout/peptide_no_h.pdb" # for visualization purposes

cat > $tmp_config_no_apbs <<- EOM
[SMIFFER]
DO_SMIF_APBS=False
EOM

cat > $tmp_config_ignore_h <<- EOM
[SMIFFER]
USE_STRUCTURE_HYDROGENS=False
DO_SMIF_APBS=False
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

rm -f $tmp_config_ignore_h $tmp_config_no_apbs
