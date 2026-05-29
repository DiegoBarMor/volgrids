#!/bin/bash
set -euo pipefail

fdata="" # PATH GOES HERE

# shellcheck disable=SC2045
for name in $(ls "$fdata"); do
    volgrids vgtools summary "$name"  | grep -E "(.*summary.*|.*min.*)" >> examples/apbs_summary/raw.txt
done
