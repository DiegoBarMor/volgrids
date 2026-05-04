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
bash tests/smiffer/run.sh
bash tests/vgtools/run.sh
bash tests/smutils/run.sh

echo "All tests completed successfully."
