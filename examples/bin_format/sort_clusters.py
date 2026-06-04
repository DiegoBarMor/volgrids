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

    cluster_ids = set(grid.arr.flatten()) - {0}
    volumes = {
        cluster_id : np.sum(grid.arr == cluster_id) for cluster_id in cluster_ids
    }.items()
    volumes = sorted(volumes, key = lambda x: x[1], reverse = True)
    cluster_ids = [pair[0] for pair in volumes if pair[1] >= VOLUME_THRESHOLD]

    new_arr = np.zeros_like(grid.arr)
    for i, cluster_id in enumerate(cluster_ids, start = 1):
        new_arr[grid.arr == cluster_id] = i

    grid.arr = new_arr
    vg.GridIO.write_bin(PATH_BIN, grid)


################################################################################
if __name__ == "__main__":
    PATH_BIN = Path(sys.argv[1])
    VOLUME_THRESHOLD = int(sys.argv[2])
    main()


################################################################################
# python3 examples/bin_format/expand_mrc2cmap.py testdata/smiffer/interfaces/prot_rna/prot.stacking.clusters.bin
