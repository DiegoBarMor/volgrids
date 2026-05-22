#!/bin/bash
set -euo pipefail

bash tests/vgtools/00_convert.sh
bash tests/vgtools/01_pack_unpack.sh
bash tests/vgtools/02_fix_cmap.sh
bash tests/vgtools/03_average.sh
bash tests/vgtools/04_summary_compare.sh
bash tests/vgtools/05_rotate_op.sh
bash tests/vgtools/06_hist.sh
