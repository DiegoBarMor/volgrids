#!/bin/bash
set -euo pipefail

folder_src=$(realpath "$(dirname "$0")")

if command -v gcc >/dev/null 2>&1; then
    CC=gcc
elif command -v clang >/dev/null 2>&1; then
    CC=clang
else
    echo "Error: neither gcc nor clang found in PATH." >&2
    exit 1
fi

echo "...>>> Compiling segmentation backend with $CC..."
"$CC" -O2 "$folder_src/main.c" -o "$folder_src/main"
