#!/bin/bash
set -eu

bash tests/env/rebuild.sh

echo
echo ">>> TEST ENV 1: Environment WITH volgrids and its dependencies installed"

root=${PWD}
fout="testdata/env"
tmp_py="smiffer.tmp.py"

cat > $tmp_py <<- EOM
import volgrids.smiffer as sm
sm.AppSmiffer.from_cli().run()
EOM

### try running it in the repo's root
python3 -W ignore $tmp_py prot testdata/_input/toy_systems/peptide.pdb -o $fout

### try running it somewhere else
cd ~
python3 -W ignore "$root"/$tmp_py prot "$root"/testdata/_input/toy_systems/peptide.pdb -o $fout

rm -f "$root"/$tmp_py
