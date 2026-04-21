#!/bin/bash
set -euo pipefail

bash volgrids/_dev/fetch_vendors.sh
python3 scripts/_update_known_configs.py volgrids
