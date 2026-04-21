import sys
from pathlib import Path

# ------------------------------------------------------------------------------
def main():
    for path in ROOT_VENDOR.rglob("*.py"):
        if path.name.startswith("_"): continue
        if path.name == "__init__.py": continue

        path.write_text(path.read_text().replace(
            f"import {NAME_VENDOR}",
            f"import volgrids._vendors.{NAME_VENDOR}"
        ))


################################################################################
if __name__ == "__main__":
    ROOT_VENDOR = Path(sys.argv[1])
    NAME_VENDOR = ROOT_VENDOR.name
    main()


################################################################################
# python3 scripts/dev/fix_vendor_imports.py volgrids/_vendors/freyacli
