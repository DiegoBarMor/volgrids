#!/bin/bash
set -eu

bash scripts/reinstall.sh

echo
echo ">>> TEST ENV 1: Environment WITH volgrids and its dependencies installed"

root=${PWD}
fout="testdata/env"

### try running it in the repo's root
volgrids smiffer prot testdata/smiffer/toy_systems/peptide.pdb -o $fout

### try running it somewhere else
cd ~
volgrids smiffer prot "$root"/testdata/smiffer/toy_systems/peptide.pdb -o $fout
