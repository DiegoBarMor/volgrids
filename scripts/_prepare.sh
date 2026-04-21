#!/bin/bash
set -euo pipefail

dir_vendors="volgrids/_vendors"

fetch_vendor_dbm() {
    local name_vendor="$1"
    dir_out="$dir_vendors/$name_vendor"

    rm -rf "$dir_out"
    git clone --depth 1 "https://github.com/DiegoBarMor/$name_vendor" $dir_vendors/_tmp
    mv "$dir_vendors/_tmp/$name_vendor" "$dir_out"
    mv $dir_vendors/_tmp/*.md "$dir_vendors/$name_vendor"/
    rm -rf $dir_vendors/_tmp

    python3 scripts/dev/fix_vendor_imports.py volgrids/_vendors/freyacli
}

mkdir -p $dir_vendors

fetch_vendor_dbm freyacli

python3 scripts/dev/update_known_configs.py volgrids
