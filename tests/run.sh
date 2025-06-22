#!/bin/bash
set -eu

tests/0-pocket_sphere.sh
tests/1-whole.sh
tests/2-traj.sh
tests/3-tools.sh
tests/4-ligand.sh

echo "All tests completed successfully."
