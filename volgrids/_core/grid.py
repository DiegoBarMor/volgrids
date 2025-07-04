import numpy as np
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Grid:
    def __init__(self, data, init_grid = True, dtype = None):
        self.ms: vg.MolecularSystem

        if isinstance(data, vg.MolecularSystem):
            self.ms = data
            self.xres, self.yres, self.zres = data.resolution
            self.xmin, self.ymin, self.zmin = data.minCoords
            self.xmax, self.ymax, self.zmax = data.maxCoords
            self.dx, self.dy, self.dz = data.deltas

        elif isinstance(data, dict):
            self.ms = None
            self.xres, self.yres, self.zres = data["resolution"]
            self.xmin, self.ymin, self.zmin = data["minCoords"]
            self.xmax, self.ymax, self.zmax = data["maxCoords"]
            self.dx, self.dy, self.dz = data["deltas"]

        if dtype is None: dtype = vg.FLOAT_DTYPE
        self.grid = np.zeros((self.xres, self.yres, self.zres), dtype = dtype) \
            if init_grid else None


    # --------------------------------------------------------------------------
    def __add__(self, other: "Grid|float|int") -> "Grid":
        grid_out = self.copy()
        # grid_out = Grid(self.ms, init_grid = False) ### [WIP] ms could be None
        if isinstance(other, Grid):
            grid_out.grid = self.grid + other.grid
            return grid_out
        try:
            grid_out.grid = self.grid + other
            return grid_out
        except TypeError:
            raise TypeError(f"Cannot add {type(other)} to Grid. Use another Grid or a numeric value.")


    # --------------------------------------------------------------------------
    def copy(self):
        vg = Grid(self.ms) ### [WIP] ms could be None
        vg.grid = np.copy(self.grid)
        return vg

    # --------------------------------------------------------------------------
    def is_empty(self):
        return np.all(self.grid == 0)

    # --------------------------------------------------------------------------
    def reshape(self, new_min: tuple[float], new_max: tuple[float], new_res: tuple[float]):
        new_xmin, new_ymin, new_zmin = new_min
        new_xmax, new_ymax, new_zmax = new_max
        new_xres, new_yres, new_zres = new_res

        self.grid = vg.Math.interpolate_3d(
            x0 = np.linspace(self.xmin, self.xmax, self.xres),
            y0 = np.linspace(self.ymin, self.ymax, self.yres),
            z0 = np.linspace(self.zmin, self.zmax, self.zres),
            data_0 = self.grid,
            new_coords = np.mgrid[
                new_xmin : new_xmax : complex(0, new_xres),
                new_ymin : new_ymax : complex(0, new_yres),
                new_zmin : new_zmax : complex(0, new_zres),
            ].T
        ).astype(vg.FLOAT_DTYPE)

        self.xmin, self.ymin, self.zmin = new_xmin, new_ymin, new_zmin
        self.xmax, self.ymax, self.zmax = new_xmax, new_ymax, new_zmax
        self.xres, self.yres, self.zres = new_xres, new_yres, new_zres
        self.dx = (self.xmax - self.xmin) / (self.xres - 1)
        self.dy = (self.ymax - self.ymin) / (self.yres - 1)
        self.dz = (self.zmax - self.zmin) / (self.zres - 1)


    # --------------------------------------------------------------------------
    @staticmethod
    def substract(grid_0: "Grid", grid_1: "Grid") -> "Grid":
        grid = grid_0.copy()
        grid.grid = grid_0.grid - grid_1.grid
        return grid

    # --------------------------------------------------------------------------
    def save_data(self, folder_out: Path, title: str):
        path_prefix = folder_out / f"{self.ms.molname}.{title}"

        if self.ms.do_traj:
            ### ignore the OUTPUT flag, CMAP is the only format that supports multiple frames
            vg.GridIO.write_cmap(f"{path_prefix}.cmap", self, f"{self.ms.molname}.{self.ms.frame:04}")
            return

        if vg.OUTPUT_FORMAT == vg.GridFormat.DX:
            vg.GridIO.write_dx(f"{path_prefix}.dx", self)
            return

        if vg.OUTPUT_FORMAT == vg.GridFormat.MRC:
            vg.GridIO.write_mrc(f"{path_prefix}.mrc", self)
            return

        if vg.OUTPUT_FORMAT == vg.GridFormat.CCP4:
            vg.GridIO.write_ccp4(f"{path_prefix}.ccp4", self)
            return

        if vg.OUTPUT_FORMAT == vg.GridFormat.CMAP:
            vg.GridIO.write_cmap(f"{path_prefix}.cmap", self, self.ms.molname)
            return

        if vg.OUTPUT_FORMAT == vg.GridFormat.CMAP_PACKED:
            vg.GridIO.write_cmap(folder_out / f"{self.ms.molname}.cmap", self, f"{self.ms.molname}.{title}")
            return

        raise ValueError(f"Unknown output format: {vg.OUTPUT_FORMAT}.")


    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
    def populate_grid(self): return


# //////////////////////////////////////////////////////////////////////////////
