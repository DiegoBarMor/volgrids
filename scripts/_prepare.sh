#!/bin/bash
set -euo pipefail

bash volgrids/_dev/fetch_vendors.sh
python3 volgrids/_dev/update_known_configs.py volgrids
