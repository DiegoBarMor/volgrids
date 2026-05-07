#!/bin/bash
set -euo pipefail

### default parameters
VERBOSE=false  # print verbose output from PDB2PQR?

help_message() {
    name_sh=$(basename "$0")
    echo "Usage: $name_sh <path_pdb> [--verbose]"
    echo "Example: $name_sh testdata/smiffer/pdb_clean/1iqj.pdb --verbose"
}
echo_red() {
    printf "\e[%sm%s\e[0m\n" "91" "$1"
}

if [[ "$#" -lt 1 ]]; then
    help_message
    exit 255
fi

path_pdb=$(realpath "$1")
shift 1

while [[ $# -gt 0 && "$1" == --* ]]; do
    case "$1" in
        --verbose)
            VERBOSE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            help_message
            exit 254
            ;;
    esac
done

if [[ ! -f "$path_pdb" ]]; then
    echo "Error: PDB file '$path_pdb' does not exist."
    exit 253
fi

folder_out=$(dirname "$path_pdb")
cd "$folder_out" || exit 252 # ------------------------------ inside output folder vvvvv

name_pdb=$(basename "$path_pdb")
path_pqr=$name_pdb.pqr
path_in="$name_pdb.in"

tmp_log=$(mktemp)
trap 'rm -f "$tmp_log"' EXIT

if [[ "$VERBOSE" == "true" ]]; then
    pdb2pqr --ff=AMBER "$name_pdb" "$path_pqr" --apbs-input "$path_in" 2>&1 | tee -a "$tmp_log" >&2
    rc=${PIPESTATUS[0]}
else
    pdb2pqr --ff=AMBER "$name_pdb" "$path_pqr" --apbs-input "$path_in" >>"$tmp_log" 2>&1
    rc=$?
fi
if [[ $rc -ne 0 ]]; then
    cat "$tmp_log" >&2
    echo_red "pdb2pqr failed (exit code $rc). See output above."
    exit "$rc"
fi

rm -f "$name_pdb.log" "$folder_out/io.mc"
