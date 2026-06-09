#!/bin/bash
set -euo pipefail

bash volgrids/_vendors/fetch_vendors.sh
python3 scripts/dev/update_known_configs.py

### delete compiled backends
rm -f volgrids/_backend/segmentation/main
