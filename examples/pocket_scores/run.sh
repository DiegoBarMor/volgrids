#!/bin/bash
set -eu

examples/pocket_scores/smifs_pocket.sh
examples/pocket_scores/smifs_whole.sh
python3 examples/pocket_scores/scores.py
