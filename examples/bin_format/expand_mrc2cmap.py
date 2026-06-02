import sys
import numpy as np
from pathlib import Path

### simulate having "volgrids" installed as a package
### this way it's not necessary to install the repo to run this script
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
### you can remove the previous two lines if volgrids is installed

import volgrids as vg

# ------------------------------------------------------------------------------
def main():
    grid = vg.GridIO.read_bin(PATH_BIN)
    grid.arr = grid.arr.astype(int)
    cluster_ids = set(grid.arr.flatten())

    for cluster_id in cluster_ids:
        if cluster_id == 0: continue

        mask = vg.Grid(grid.box, init_grid = False)
        mask.arr = grid.arr == cluster_id

        volume = np.sum(mask.arr)
        if volume < VOLUME_THRESHOLD: continue

        print(f"Volume cluster {cluster_id}: {volume}")
        vg.GridIO.write_cmap(PATH_CMAP, mask, key = f"cluster.{cluster_id}")


################################################################################
if __name__ == "__main__":
    PATH_BIN  = Path(sys.argv[1])
    PATH_CMAP = PATH_BIN.with_suffix(".cmap")
    VOLUME_THRESHOLD = 200
    main()


################################################################################
# python3 examples/bin_format/expand_mrc2cmap.py testdata/smiffer/interfaces/prot_rna/prot.stacking.clusters.bin
