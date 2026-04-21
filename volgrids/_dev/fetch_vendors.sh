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

    echo "Fetched $name_vendor $(cat "$dir_vendors/$name_vendor/_version.py")"

    python3 - "$dir_vendors/$name_vendor" <<'PYCODE'
import sys
from pathlib import Path

root_vendor = Path(sys.argv[1])
name_vendor = root_vendor.name
for path in root_vendor.rglob("*.py"):
    if path.name.startswith("_"): continue
    if path.name == "__init__.py": continue
    path.write_text(path.read_text().replace(
        f"import {name_vendor}",
        f"import volgrids._vendors.{name_vendor}"
    ))
PYCODE
}

mkdir -p $dir_vendors

fetch_vendor_dbm freyacli
