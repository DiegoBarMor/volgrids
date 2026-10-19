#!/bin/bash
set -euo pipefail

echo
echo ">>> TEST SMUTILS 0: RNA non-interacting residues"

fpdb="testdata/smiffer/pdb_clean"

_assert_result() {
    local smu_output; read -r smu_output;
    local expected=$1
    local command=$2
    if [ "$smu_output" != "$expected" ]; then
        echo "XXX Failed test for '$command'"
        echo "Expected: $expected"
        echo "Output: $smu_output $(</dev/stdin)"
        exit 1
    fi
}

run_test() {
    local name_pdb=$1
    local command=$2
    local expected=$3

    python3 volgrids smutils "$command" "$fpdb/$name_pdb.pdb" | _assert_result "$expected" "$command"
    echo ">>> Passed '$command' test for '$name_pdb'"
}

run_test 1akx res_nobp  "A.23 A.25 A.30 A.31 A.32 A.33 A.34 A.35"
run_test 1akx res_nostk "A.16 A.19 A.22 A.23 A.25 A.26 A.31 A.32 A.33 A.34 A.35 A.46"

run_test 1i9v res_nobp  "A.14 A.15 A.16 A.17 A.18 A.20 A.21 A.26 A.31 A.32 A.33 A.34 A.35 A.36 A.38 A.39 A.4 A.44 A.45 A.46 A.47 A.48 A.55 A.57 A.58 A.59 A.60 A.69 A.73 A.74 A.75 A.76 A.8 A.9"
run_test 1i9v res_nostk "A.1 A.10 A.15 A.16 A.17 A.19 A.20 A.33 A.34 A.36 A.38 A.39 A.40 A.45 A.47 A.49 A.50 A.55 A.56 A.60 A.61 A.75 A.76 A.8 A.9"

run_test 2esj res_nobp  "A.16 A.17 A.19 A.5 A.7 B.28 B.30 B.39 B.40 B.42"
run_test 2esj res_nostk "A.14 A.15 A.16 A.17 A.2 A.20 A.22 B.25 B.26 B.37 B.38 B.39 B.40 B.45"

run_test 4f8u res_nobp  "A.17 A.18 A.2 A.20 A.6 A.8 B.29 B.31 B.40 B.41 B.43"
run_test 4f8u res_nostk "A.12 A.17 A.18 A.2 A.23 A.4 B.26 B.30 B.40 B.41 B.44 B.46"

run_test 5bjo res_nobp  "E.10 E.11 E.12 E.13 E.14 E.15 E.16 E.19 E.21 E.22 E.23 E.24 E.25 E.26 E.27 E.29 E.30 E.31 E.6 E.7 E.8 Y.10 Y.11 Y.12 Y.13 Y.14 Y.15 Y.16 Y.19 Y.21 Y.22 Y.23 Y.24 Y.25 Y.26 Y.27 Y.29 Y.30 Y.31 Y.6 Y.7 Y.8"
run_test 5bjo res_nostk "E.1 E.11 E.12 E.14 E.15 E.16 E.18 E.19 E.20 E.24 E.25 E.28 E.31 E.36 E.8 Y.1 Y.10 Y.11 Y.12 Y.13 Y.14 Y.16 Y.18 Y.19 Y.20 Y.22 Y.24 Y.25 Y.28 Y.3 Y.31 Y.36 Y.8"

run_test 5kx9 res_nobp  "X.10 X.11 X.12 X.13 X.14 X.18 X.19 X.2 X.20 X.21 X.22 X.23 X.24 X.28 X.37 X.38 X.39 X.40 X.47 X.48 X.53 X.9 Y.100 Y.101 Y.102 Y.103 Y.104 Y.111 Y.61 Y.62 Y.63 Y.69 Y.70 Y.71 Y.72 Y.73 Y.74 Y.75 Y.85 Y.89 Y.90 Y.91 Y.92 Y.98 Y.99"
run_test 5kx9 res_nostk "X.1 X.10 X.12 X.13 X.14 X.19 X.20 X.22 X.23 X.24 X.38 X.39 X.46 X.48 X.49 X.52 X.53 Y.100 Y.104 Y.112 Y.57 Y.61 Y.62 Y.63 Y.70 Y.71 Y.72 Y.74 Y.85 Y.89 Y.90 Y.91 Y.94"

run_test 6tf3 res_nobp  "A.10 A.11 A.12 A.15 A.20 A.24 A.28 A.29 A.30 A.31 A.35 A.41 A.48 A.5 A.8 A.9"
run_test 6tf3 res_nostk "A.1 A.10 A.11 A.12 A.28 A.29 A.35 A.36 A.37 A.38 A.52 A.8 A.9"

run_test 7oax0 res_nobp  "A.11 A.12 A.13 A.14 A.17 A.22 A.23 A.24 A.25 A.30 A.31 A.32 A.33 A.34 A.35 A.36 A.37 A.38 A.41 A.42 A.43 A.45 A.52 A.8 A.9 B.11 B.12 B.13 B.14 B.17 B.22 B.23 B.24 B.25 B.30 B.31 B.32 B.33 B.34 B.35 B.36 B.37 B.38 B.41 B.42 B.43 B.45 B.52 B.8 B.9 C.11 C.12 C.13 C.14 C.17 C.22 C.23 C.24 C.25 C.30 C.31 C.32 C.33 C.34 C.35 C.36 C.37 C.38 C.41 C.42 C.43 C.45 C.52 C.8 C.9 D.11 D.12 D.13 D.14 D.17 D.22 D.23 D.24 D.25 D.30 D.31 D.32 D.33 D.34 D.35 D.36 D.37 D.38 D.41 D.42 D.43 D.45 D.52 D.8 D.9"
run_test 7oax0 res_nostk "A.10 A.11 A.13 A.14 A.2 A.23 A.24 A.25 A.31 A.32 A.33 A.34 A.35 A.37 A.38 A.42 A.43 A.52 A.8 A.9 B.10 B.11 B.13 B.14 B.2 B.23 B.24 B.25 B.31 B.32 B.33 B.34 B.35 B.36 B.37 B.38 B.43 B.52 B.8 B.9 C.10 C.11 C.13 C.14 C.2 C.23 C.24 C.25 C.27 C.31 C.32 C.33 C.34 C.35 C.37 C.38 C.42 C.43 C.52 C.8 C.9 D.10 D.11 D.13 D.14 D.2 D.23 D.24 D.25 D.31 D.32 D.33 D.34 D.35 D.38 D.43 D.52 D.8 D.9"

run_test 7oax1 res_nobp  "A.11 A.12 A.13 A.14 A.17 A.22 A.23 A.24 A.25 A.30 A.31 A.32 A.33 A.34 A.35 A.36 A.37 A.38 A.41 A.42 A.43 A.45 A.52 A.8 A.9 B.11 B.12 B.13 B.14 B.17 B.22 B.23 B.24 B.25 B.30 B.31 B.32 B.33 B.34 B.35 B.36 B.37 B.38 B.41 B.42 B.43 B.45 B.52 B.8 B.9 C.11 C.12 C.13 C.14 C.17 C.22 C.23 C.24 C.25 C.30 C.31 C.32 C.33 C.34 C.35 C.36 C.37 C.38 C.41 C.42 C.43 C.45 C.52 C.8 C.9 D.11 D.12 D.13 D.14 D.17 D.22 D.23 D.24 D.25 D.30 D.31 D.32 D.33 D.34 D.35 D.36 D.37 D.38 D.41 D.42 D.43 D.45 D.52 D.8 D.9"
run_test 7oax1 res_nostk "A.10 A.11 A.13 A.14 A.2 A.23 A.24 A.25 A.31 A.32 A.33 A.34 A.35 A.37 A.38 A.42 A.43 A.52 A.8 A.9 B.10 B.11 B.13 B.14 B.2 B.23 B.24 B.25 B.31 B.32 B.33 B.34 B.35 B.36 B.37 B.38 B.43 B.52 B.8 B.9 C.10 C.11 C.13 C.14 C.2 C.23 C.24 C.25 C.27 C.31 C.32 C.33 C.34 C.35 C.37 C.38 C.42 C.43 C.52 C.8 C.9 D.10 D.11 D.13 D.14 D.2 D.23 D.24 D.25 D.31 D.32 D.33 D.34 D.35 D.38 D.43 D.52 D.8 D.9"

run_test 8eyv res_nobp  "A.12 A.17 A.18 A.19 A.20 A.21 A.22 A.23 A.24 A.25 A.27 A.28 A.29 A.30 A.31 A.32 A.33 A.34 A.35 A.36 A.42 B.12 B.17 B.18 B.19 B.20 B.21 B.22 B.23 B.24 B.25 B.27 B.28 B.29 B.30 B.31 B.32 B.33 B.34 B.35 B.36 B.42"
run_test 8eyv res_nostk "A.1 A.16 A.17 A.18 A.19 A.2 A.23 A.26 A.27 A.28 A.29 A.31 A.34 A.36 A.37 A.38 A.49 B.1 B.11 B.12 B.16 B.17 B.18 B.19 B.20 B.22 B.23 B.24 B.26 B.27 B.28 B.31 B.33 B.34 B.36 B.37 B.38 B.47 B.49"
