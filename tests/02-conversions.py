from pathlib import Path
import volpot as vp

FOLDER_TEST = Path("data/tests/02-conversions")
FOLDER_TEST.mkdir(parents = True, exist_ok = True)

##### Input files
PATH_MRC_ORIG = FOLDER_TEST / "1iqj.stk.mrc"
PATH_DX_ORIG  = FOLDER_TEST / "1iqj.stk.dx"

##### Output files
PATH_MRC_TO_MRC = FOLDER_TEST / "mrc2mrc.mrc"
PATH_MRC_TO_DX  = FOLDER_TEST / "mrc2dx.dx"
PATH_DX_TO_MRC  = FOLDER_TEST / "dx2mrc.mrc"
PATH_DX_TO_DX   = FOLDER_TEST / "dx2dx.dx"

##### Test operations
print("\n>>> TEST 02: Conversions between MRC and DX formats")

print(f"...>>> Converting {PATH_MRC_ORIG}")
grid_mrc = vp.read_mrc(PATH_MRC_ORIG)
vp.write_mrc(PATH_MRC_TO_MRC, grid_mrc)
vp.write_dx(PATH_MRC_TO_DX, grid_mrc)

print(f"...>>> Converting {PATH_DX_ORIG}")
grid_dx = vp.read_dx(PATH_DX_ORIG)
vp.write_mrc(PATH_DX_TO_MRC, grid_dx)
vp.write_dx(PATH_DX_TO_DX, grid_dx)

# vmd data/pdb/1iqj.pdb data/tests/1iqj.dx2dx.dx data/tests/1iqj.mrc2dx.dx data/tests/1iqj.stk.dx data/tests/1iqj.dx2mrc.mrc data/tests/1iqj.mrc2mrc.mrc data/tests/1iqj.stk.mrc
