#!/bin/bash
set -euo pipefail

bash tests/smiffer/00_toy_systems.sh
bash tests/smiffer/01_pocket_sphere.sh
bash tests/smiffer/02_whole.sh
bash tests/smiffer/03_rna_hbonds.sh
bash tests/smiffer/04_ligand.sh
bash tests/smiffer/05_traj.sh
bash tests/smiffer/06_cavities.sh
# bash tests/smiffer/07_box_csv.sh
