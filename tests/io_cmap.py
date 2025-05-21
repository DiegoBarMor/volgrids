from pathlib import Path
import volpot as vp

FOLDER_TEST = Path("data/potentials-test")
FOLDER_OUT = FOLDER_TEST / "cmap"
FOLDER_OUT.mkdir(parents = True, exist_ok = True)

##### Paths Unpacking...
PATH_CMAP_IN = FOLDER_TEST / "cell15_timeseries.cmap"
PATH_MRC_OUT_UNPACKED_PREFFIX = FOLDER_OUT / "cell15"

##### Paths Packing...
PATHS_MRC_IN = [
    FOLDER_TEST / "1iqj.hba.mrc",
    FOLDER_TEST / "1iqj.hbd.mrc",
    FOLDER_TEST / "1iqj.phi.mrc",
    FOLDER_TEST / "1iqj.pho.mrc",
    FOLDER_TEST / "1iqj.stk.mrc",
]
PATH_CMAP_PACKED = FOLDER_OUT / "1iqj.cmap"

# ------------------------------------------------------------------------------
for key in vp.get_cmap_keys(PATH_CMAP_IN):
    print(">>> Unpacking cmap key:", key)
    grid = vp.read_cmap(PATH_CMAP_IN, key)
    vp.write_mrc(f"{PATH_MRC_OUT_UNPACKED_PREFFIX}.{key}.mrc", grid)

for path_mrc in PATHS_MRC_IN:
    print(">>> Packing mrc:", path_mrc)
    grid = vp.read_mrc(path_mrc)
    vp.write_cmap(PATH_CMAP_PACKED, grid, path_mrc.stem)
