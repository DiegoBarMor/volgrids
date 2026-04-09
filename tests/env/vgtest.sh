#!/bin/bash
set -eu



echo
echo ">>> TEST ENV 0: Environment WITHOUT volgrids installed, just its dependencies"

root=${PWD}
fout="testdata/env"

### try running it in the repo's root
python3 volgrids smiffer prot testdata/smiffer/toy_systems/peptide.pdb -o $fout

### try running it somewhere else
cd ~
python3 "$root"/volgrids smiffer prot "$root"/testdata/smiffer/toy_systems/peptide.pdb -o $fout
