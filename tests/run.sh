#!/bin/bash
set -eu

tests/smiffer/toy_systems.sh
tests/smiffer/pocket_sphere.sh
tests/smiffer/whole.sh
tests/smiffer/traj.sh
tests/smiffer/ligand.sh

tests/vgtools/tools.sh

echo "All tests completed successfully."
