#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 0: Toy systems"

fout="testdata/smiffer/toy_systems"

conf_no_apbs="SMIF_APBS=False"
conf_ignore_h="$conf_no_apbs SMIF_USE_HYDROGENS=False"
conf_equilateral="$conf_no_apbs BOX_FORCE_EQUILATERAL=True"

python3 volgrids smiffer $fout/peptide_no_h.pdb    -o $fout -c "$conf_ignore_h"
python3 volgrids smiffer $fout/peptide.pdb         -o $fout -c $conf_no_apbs
python3 volgrids smiffer $fout/guanine.pdb         -o $fout -c $conf_no_apbs
python3 volgrids smiffer $fout/ribose_gua_no_h.pdb -o $fout -c "$conf_ignore_h"
python3 volgrids smiffer $fout/ribose_gua.pdb      -o $fout -c $conf_no_apbs
python3 volgrids smiffer $fout/uuu.pdb             -o $fout -c $conf_no_apbs

python3 volgrids smiffer $fout/all_arg.pdb -o $fout -c $conf_no_apbs
python3 volgrids smiffer $fout/all_asn.pdb -o $fout -c $conf_no_apbs
python3 volgrids smiffer $fout/all_cyt.pdb -o $fout -c $conf_no_apbs
python3 volgrids smiffer $fout/all_ump.pdb -o $fout -c $conf_no_apbs

cp $fout/all_cyt.pdb $fout/all_cyt_equilateral.pdb
# shellcheck disable=SC2086
python3 volgrids smiffer $fout/all_cyt_equilateral.pdb -o $fout -c $conf_equilateral # tests also expansion of -c arguments
rm -f $fout/all_cyt_equilateral.pdb
