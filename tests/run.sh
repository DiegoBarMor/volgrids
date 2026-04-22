#!/bin/bash
set -eu

##### Generate input data (if needed)
if [[ ! -d testdata/smiffer/apbs ]]; then
    echo "Generating input data..."
    bash tests/_gen_input.sh
else
    echo "Input data already exists, skipping generation."
fi


##### Run tests
bash tests/env/vgtest.sh

bash tests/smiffer/00_toy_systems.sh
bash tests/smiffer/01_pocket_sphere.sh
bash tests/smiffer/02_whole.sh
bash tests/smiffer/03_rna_hbonds.sh
bash tests/smiffer/04_ligand.sh
bash tests/smiffer/05_traj.sh
bash tests/smiffer/06_cavities.sh

bash tests/vgtools/00_convert.sh
bash tests/vgtools/01_pack_unpack.sh
bash tests/vgtools/02_fix_cmap.sh
bash tests/vgtools/03_average.sh
bash tests/vgtools/04_summary_compare.sh
bash tests/vgtools/05_rotate_op.sh

echo "All tests completed successfully."
