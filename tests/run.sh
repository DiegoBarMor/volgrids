#!/bin/bash
set -eu

tests/env/vgtest.sh

tests/smiffer/toy_systems.sh
tests/smiffer/pocket_sphere.sh
tests/smiffer/whole.sh
tests/smiffer/traj.sh
tests/smiffer/ligand.sh

tests/vgtools/convert.sh
tests/vgtools/pack_unpack.sh
tests/vgtools/fix_cmap.sh
tests/vgtools/compare.sh

echo "All tests completed successfully."
