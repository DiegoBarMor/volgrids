#!/bin/bash
set -euo pipefail

### Automatic calculation of electrostatic potential grids using APBS.
### The grid will be saved in the same folder as the input PDB file,
### with the same name but with extension .dx or .mrc (if --mrc is used).
### All intermediate files will be removed (PQR can be conserved by using --pqr).
###
### Requires: pdb2pqr, apbs
### Optional: python3 with volgrids installed (for MRC conversion)
### Usage: apbs.sh <path_pdb> [--mrc] [--pqr] [--verbose]
### Example: apbs.sh testdata/smiffer/pdb-nosolv/1iqj.pdb --mrc --pqr --verbose


### default parameters
CONVERT_TO_MRC=false # convert the output to MRC format instead of DX?
KEEP_PQR=false       # keep the intermediate PQR file?
VERBOSE=false        # print verbose output from APBS?

help_message() {
    name_sh=$(basename "$0")
    echo "Usage: $name_sh <path_pdb> [--mrc] [--verbose]"
    echo "Example: $name_sh testdata/smiffer/pdb-nosolv/1iqj.pdb --mrc --verbose"
}

if [[ "$#" -lt 1 ]]; then
    help_message
    exit 1
fi

path_pdb=$(realpath "$1")
shift 1

while [[ $# -gt 0 && "$1" == --* ]]; do
    case "$1" in
        --mrc)
            CONVERT_TO_MRC=true
            shift
            ;;
        --pqr)
            KEEP_PQR=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            help_message
            exit 1
            ;;
    esac
done

cwd=$(pwd)

if [[ ! -f "$path_pdb" ]]; then
    echo "Error: PDB file '$path_pdb' does not exist."
    exit 1
fi

folder_out=$(dirname "$path_pdb")
cd "$folder_out" # ------------------------------ inside output folder vvvvv

name_pdb=$(basename "$path_pdb")
path_pqr=$name_pdb.pqr
path_in="$name_pdb.in"

tmp_log=$(mktemp)
trap 'rm -f "$tmp_log"' EXIT

### Create APBS input file
if [[ "$VERBOSE" == "true" ]]; then
    pdb2pqr --ff=AMBER "$name_pdb" "$path_pqr" --apbs-input "$path_in" 2>&1 | tee -a "$tmp_log" >&2
    rc=${PIPESTATUS[0]}
else
    pdb2pqr --ff=AMBER "$name_pdb" "$path_pqr" --apbs-input "$path_in" >>"$tmp_log" 2>&1
    rc=$?
fi
if [[ $rc -ne 0 ]]; then
    echo "pdb2pqr failed (exit code $rc). Full output in: $tmp_log" >&2
    cat "$tmp_log" >&2
    exit "$rc"
fi

### Run APBS
if [[ "$VERBOSE" == "true" ]]; then
    apbs "$path_in" 2>&1 | tee "$tmp_log" >&2
    rc=${PIPESTATUS[0]}
else
    apbs "$path_in" >"$tmp_log" 2>&1
    rc=$?
fi
if [[ $rc -ne 0 ]]; then
    echo "APBS failed (exit code $rc). Full output in: $tmp_log" >&2
    cat "$tmp_log" >&2
    exit "$rc"
fi

log=$(grep -m1 "Writing potential to" "$tmp_log" || true)
if [[ -z "$log" ]]; then
    echo "Could not find 'Writing potential to' in APBS output. Full output below:" >&2
    cat "$tmp_log" >&2
    exit 1
fi
path_grid=$(echo "$log" | awk '{print $NF}')

mv "$path_grid" "$name_pdb.dx"

rm -f "$name_pdb.in" "$name_pdb.log" "$folder_out/io.mc"
if [[ "$KEEP_PQR" == "false" ]]; then
    rm -f "$path_pqr"
fi

cd "$cwd"  # ------------------------------------ back to original folder vvvvv

if [[ "$CONVERT_TO_MRC" == "true" ]]; then
    root_volgrids=$(dirname "$(dirname "$(realpath "$0")")")
    ### if this sh script is moved somewhere else and the volgrids package is installed
    ### the command below can be simplified to just "volgrids vgtools convert ..."
    python3 "$root_volgrids" vgtools convert "$folder_out/$name_pdb.dx" --mrc "$folder_out/$name_pdb.mrc"
    rm -f "$folder_out/$name_pdb.dx"
fi
