#!/bin/bash
set -eu

##### Fetch testdata (if needed)
if [[ ! -d testdata/smiffer/apbs ]]; then
    git clone --depth 1 https://github.com/DiegoBarMor/volgrids-testdata testdata
fi

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
