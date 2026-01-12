#!/bin/bash
set -eu

echo
echo ">>> TEST ENV 0: Environment WITHOUT volgrids installed, just its dependencies"

root=${PWD}
fout="testdata/env"

### try running it in the repo's root
python3 smiffer.py prot testdata/_raw_input/toy_systems/peptide.pdb -o $fout

### try running it somewhere else
cd ~
python3 "$root"/smiffer.py prot "$root"/testdata/_raw_input/toy_systems/peptide.pdb -o $fout

cd "$root"
rm -rf $fout
