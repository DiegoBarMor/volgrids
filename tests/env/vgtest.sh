#!/bin/bash
set -eu

echo
echo ">>> TEST ENV 0: Environment WITHOUT volgrids installed, just its dependencies"

root=${PWD}
fout="testdata/env"

### try running it in the repo's root
python3 run/smiffer.py prot testdata/_input/toy_systems/peptide.pdb -o $fout

### try running it somewhere else
cd ~
python3 "$root"/run/smiffer.py prot "$root"/testdata/_input/toy_systems/peptide.pdb -o $fout

cd "$root"
rm -rf $fout
