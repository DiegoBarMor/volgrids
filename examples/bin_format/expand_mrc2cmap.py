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
    grid = vg.Grid.load(PATH_BIN, vg.GridFormat.BIN)
    grid.arr = grid.arr.astype(int)
    cluster_ids = set(grid.arr.flatten()) - {0}

    for cluster_id in cluster_ids:
        mask = vg.Grid(grid.box, init_grid = False)
        mask.arr = grid.arr == cluster_id

        volume = np.sum(mask.arr)

        print(f"Volume cluster {cluster_id}: {volume}")
        mask.save(PATH_CMAP, vg.GridFormat.CMAP, key = f"cluster.{cluster_id:04}")


################################################################################
if __name__ == "__main__":
    PATH_BIN  = Path(sys.argv[1])
    PATH_CMAP = PATH_BIN.with_suffix(".cmap")
    main()


################################################################################
# python3 examples/bin_format/expand_mrc2cmap.py testdata/smiffer/interfaces/prot_rna/prot.stk.clusters.bin
