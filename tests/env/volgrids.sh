#!/bin/bash
set -eu

bash tests/env/rebuild.sh

echo
echo ">>> TEST ENV 1: Environment WITH volgrids and its dependencies installed"

root=${PWD}
fout="testdata/env"

### try running it in the repo's root
volgrids smiffer prot testdata/_raw_input/toy_systems/peptide.pdb -o $fout

### try running it somewhere else
cd ~
volgrids smiffer prot "$root"/testdata/_raw_input/toy_systems/peptide.pdb -o $fout

