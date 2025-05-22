from pathlib import Path
import volpot as vp

FOLDER_TEST = Path("data/tests/03-cmap")
FOLDER_TEST.mkdir(parents = True, exist_ok = True)

##### Input files
PATH_CMAP_IN = FOLDER_TEST / "1iqj.cmap"
PATHS_MRC_IN = [
    FOLDER_TEST / "2esj.hba.mrc",
    FOLDER_TEST / "2esj.hbd.mrc",
    FOLDER_TEST / "2esj.phi.mrc",
    FOLDER_TEST / "2esj.pho.mrc",
    FOLDER_TEST / "2esj.stk.mrc",
]

##### Output files
PATH_CMAP_PACKED = FOLDER_TEST / "2esj.cmap"

##### Test operations
print("\n>>> TEST 03: CMAP packing / unpacking")

for key in vp.get_cmap_keys(PATH_CMAP_IN):
    print("...>>> Unpacking cmap key:", key)
    grid = vp.read_cmap(PATH_CMAP_IN, key)
    vp.write_mrc(FOLDER_TEST / f"{key}.mrc", grid)

for path_mrc in PATHS_MRC_IN:
    print("...>>> Packing mrc:", path_mrc)
    grid = vp.read_mrc(path_mrc)
    vp.write_cmap(PATH_CMAP_PACKED, grid, path_mrc.stem)
