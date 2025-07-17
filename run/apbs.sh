#!/bin/bash
set -eu

if [[ "$#" -ne 2 ]]; then
    echo "Usage: $0 <path_pdb> <folder_output>"
    echo "Example: $0 testdata/_input/pdb-nosolv/1iqj.pdb testdata/_input/apbs"
    exit 1
fi

### set this to "true" if the output should be in MRC format instead of DX
### MRC conversion is provided by the `vgtools` converter from this repo
CONVERT_TO_MRC=true


path_pdb="$1"
folder_out="$2"
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

pdb2pqr --ff=AMBER "$name_pdb" "$path_pqr" --apbs-input "$path_in"
log=$(apbs "$path_in" | tee /dev/stderr | grep "  Writing potential to ")
path_grid=$(echo "$log" | awk '{print $NF}')

mv "$path_grid" "$name_pdb.dx"

cd "$cwd"  # ------------------------------------ back to root folder vvvvv

preffix="$folder_out/$name_pdb"
rm -f "$preffix" "$preffix.in" "$preffix.log" "$preffix.pqr" "$folder_out/io.mc"
if [[ "$CONVERT_TO_MRC" == "true" ]]; then
    python3 run/vgtools.py convert "$preffix.dx" --mrc "$preffix.mrc"
    rm -f "$preffix.dx"
fi
