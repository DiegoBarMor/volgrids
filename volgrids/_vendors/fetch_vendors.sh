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
}
fix_vendor_imports() {
    local name_vendor="$1"
    local name_imports="$2"
    python3 - "$dir_vendors/$name_vendor" <<PYCODE
import sys
from pathlib import Path

root_vendor = Path(sys.argv[1])
for path in root_vendor.rglob("*.py"):
    if path.name.startswith("_"): continue
    if path.name == "__init__.py": continue
    path.write_text(path.read_text().replace(
        "import $name_imports",
        "import volgrids._vendors.$name_imports"
    ))
PYCODE
}

mkdir -p $dir_vendors

fetch_vendor_dbm freyacli
fetch_vendor_dbm molutils

fix_vendor_imports freyacli freyacli
fix_vendor_imports molutils molutils
fix_vendor_imports molutils freyacli
