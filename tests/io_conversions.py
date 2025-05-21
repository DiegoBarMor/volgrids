from pathlib import Path
import volpot as vp

FOLDER_TEST = Path("data/potentials-test")

PATH_MRC_ORIG = FOLDER_TEST / "1iqj.stk.mrc"
PATH_DX_ORIG  = FOLDER_TEST / "1iqj.stk.dx"

PATH_MRC_TO_MRC = FOLDER_TEST / "1iqj.mrc2mrc.mrc"
PATH_MRC_TO_DX  = FOLDER_TEST / "1iqj.mrc2dx.dx"
PATH_DX_TO_MRC  = FOLDER_TEST / "1iqj.dx2mrc.mrc"
PATH_DX_TO_DX   = FOLDER_TEST / "1iqj.dx2dx.dx"

grid_mrc = vp.Grid.read_mrc(PATH_MRC_ORIG)
vp.write_mrc(PATH_MRC_TO_MRC, grid_mrc)
vp.write_dx(PATH_MRC_TO_DX, grid_mrc)

grid_dx = vp.Grid.read_dx(PATH_DX_ORIG)
vp.write_mrc(PATH_DX_TO_MRC, grid_dx)
vp.write_dx(PATH_DX_TO_DX, grid_dx)

# vmd data/pdb/1iqj.pdb data/potentials-test/1iqj.dx2dx.dx data/potentials-test/1iqj.mrc2dx.dx data/potentials-test/1iqj.stk.dx data/potentials-test/1iqj.dx2mrc.mrc data/potentials-test/1iqj.mrc2mrc.mrc data/potentials-test/1iqj.stk.mrc
