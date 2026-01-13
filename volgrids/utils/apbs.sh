#!/bin/bash
set -euo pipefail

### convert the output to MRC format instead of DX
### MRC conversion is provided by the `vgtools` converter from this repo
CONVERT_TO_MRC=false
VERBOSE=false

if [[ "$#" -lt 3 ]]; then
    echo "Usage: $0 <path_pdb> <folder_output> [--mrc] [--verbose]"
    echo "Example: $0 testdata/smiffer/pdb-nosolv/1iqj.pdb testdata/smiffer/apbs --mrc --verbose"
    exit 1
fi

path_pdb="$1"
folder_out="$2"
shift 2

while [[ $# -gt 0 && "$1" == --* ]]; do
    case "$1" in
        --mrc)
            CONVERT_TO_MRC=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 <path_pdb> <folder_output> [--mrc] [--verbose]"
            exit 1
            ;;
    esac
done

cwd=$(pwd)

if [[ ! -f "$path_pdb" ]]; then
    echo "Error: PDB file '$path_pdb' does not exist."
    exit 1
fi
if [[ ! -d "$folder_out" ]]; then
    echo "Error: Output folder '$folder_out' does not exist."
    exit 1
fi

cp "$path_pdb" "$folder_out"

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

cd "$cwd"  # ------------------------------------ back to original folder vvvvv

preffix="$folder_out/$name_pdb"
rm -f "$preffix" "$preffix.in" "$preffix.log" "$preffix.pqr" "$folder_out/io.mc"

if [[ "$CONVERT_TO_MRC" == "true" ]]; then
    python3 vgtools.py convert "$preffix.dx" --mrc "$preffix.mrc"
    rm -f "$preffix.dx"
fi
