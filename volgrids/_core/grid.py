import numpy as np
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Grid:
    def __init__(self, ms: "vg.MolSystem", init_grid = True, dtype = None):
        self.ms = ms
        self.dtype: type = vg.FLOAT_DTYPE if dtype is None else dtype
        self.arr: np.ndarray|None = np.zeros(ms.resolution, dtype = self.dtype) if init_grid else None
        self.fmt: vg.GridFormat = None


    # --------------------------------------------------------------------------
    def __add__(self, other: "Grid|float|int") -> "Grid":
        obj = Grid(self.ms, init_grid = False)
        if isinstance(other, Grid):
            obj.arr = self.arr + other.arr
            return obj
        try:
            obj.arr = self.arr + other
            return obj
        except TypeError:
            raise TypeError(f"Cannot add {type(other)} to Grid. Use another Grid or a numeric value.")


    # --------------------------------------------------------------------------
    def __sub__(self, other: "Grid|float|int") -> "Grid":
        obj = Grid(self.ms, init_grid = False)
        if isinstance(other, Grid):
            obj.arr = self.arr - other.arr
            return obj
        try:
            obj.arr = self.arr - other
            return obj
        except TypeError:
            raise TypeError(f"Cannot substract {type(other)} from Grid. Use another Grid or a numeric value.")


    # --------------------------------------------------------------------------
    def __abs__(self) -> "Grid":
        obj = Grid(self.ms, init_grid = False)
        obj.arr = np.abs(self.arr)
        return obj


    # --------------------------------------------------------------------------
    @classmethod
    def reverse(cls, other: "Grid") -> "Grid":
        """Return a new Grid with the reversed values of the other Grid.
        For boolean grids, the reverse is the logical not.
        For numeric grids, the reverse is the negation of the values.
        """
        obj = cls(other.ms, init_grid = False)
        obj.arr = np.logical_not(other.arr) if (other.dtype == bool) else -other.arr
        return obj

    # -------------------------------------------------------------------------- GETTERS
    def xres(self): return self.ms.resolution[0]
    def yres(self): return self.ms.resolution[1]
    def zres(self): return self.ms.resolution[2]
    def xmin(self): return self.ms.min_coords[0]
    def ymin(self): return self.ms.min_coords[1]
    def zmin(self): return self.ms.min_coords[2]
    def xmax(self): return self.ms.max_coords[0]
    def ymax(self): return self.ms.max_coords[1]
    def zmax(self): return self.ms.max_coords[2]
    def   dx(self): return self.ms.deltas[0]
    def   dy(self): return self.ms.deltas[1]
    def   dz(self): return self.ms.deltas[2]

    def npoints(self): return self.xres() * self.yres() * self.zres()

    # --------------------------------------------------------------------------
    def copy(self):
        obj = Grid(self.ms, init_grid = False)
        obj.arr = np.copy(self.arr)
        return obj


    # --------------------------------------------------------------------------
    def is_empty(self):
        return np.all(self.arr == 0)


    # --------------------------------------------------------------------------
    def reshape(self, new_min: tuple[float], new_max: tuple[float], new_res: tuple[float]):
        new_xmin, new_ymin, new_zmin = new_min
        new_xmax, new_ymax, new_zmax = new_max
        new_xres, new_yres, new_zres = new_res

        self.arr = vg.Math.interpolate_3d(
            x0 = np.linspace(self.xmin(), self.xmax(), self.xres()),
            y0 = np.linspace(self.ymin(), self.ymax(), self.yres()),
            z0 = np.linspace(self.zmin(), self.zmax(), self.zres()),
            data_0 = self.arr,
            new_coords = np.mgrid[
                new_xmin : new_xmax : complex(0, new_xres),
                new_ymin : new_ymax : complex(0, new_yres),
                new_zmin : new_zmax : complex(0, new_zres),
            ].T
        ).astype(vg.FLOAT_DTYPE)

        self.ms.min_coords = np.array([new_xmin, new_ymin, new_zmin])
        self.ms.max_coords = np.array([new_xmax, new_ymax, new_zmax])
        self.ms.resolution = np.array([new_xres, new_yres, new_zres])
        self.ms.deltas = (self.ms.max_coords - self.ms.min_coords) / (self.ms.resolution - 1)


    # --------------------------------------------------------------------------
    def has_equivalent_box(self, other: "Grid") -> bool:
        return not np.any(
            (self.ms.resolution - other.ms.resolution) +\
            (self.ms.min_coords - other.ms.min_coords) +\
            (self.ms.max_coords - other.ms.max_coords)
        )


    # --------------------------------------------------------------------------
    def reshape_as(self, other: "Grid"):
        self.reshape(
            new_min = other.ms.min_coords,
            new_max = other.ms.max_coords,
            new_res = other.ms.resolution,
        )


    # --------------------------------------------------------------------------
    def save_data(self, folder_out: Path, title: str):
        def clear_cmap(path_out: Path) -> Path:
            if vg.REMOVE_OLD_CMAP_OUTPUT:
                path_out.unlink(missing_ok = True)
            return path_out

        path_prefix = folder_out / f"{self.ms.molname}.{title}"

        if self.ms.do_traj: # ignores the GRID_FORMAT_OUTPUT config -> CMAP is the only format that supports multiple frames
            path_cmap = clear_cmap(folder_out / f"{self.ms.molname}.{title}.cmap")
            vg.REMOVE_OLD_CMAP_OUTPUT = False
            vg.GridIO.write_cmap(path_cmap, self, f"{self.ms.molname}.{self.ms.frame:04}")
            return

        gf = vg.GridFormat.from_str(vg.GRID_FORMAT_OUTPUT)

        if gf == vg.GridFormat.DX:
            vg.GridIO.write_dx(f"{path_prefix}.dx", self)
            return

        if gf == vg.GridFormat.MRC:
            vg.GridIO.write_mrc(f"{path_prefix}.mrc", self)
            return

        if gf == vg.GridFormat.CCP4:
            vg.GridIO.write_ccp4(f"{path_prefix}.ccp4", self)
            return

        if gf == vg.GridFormat.CMAP:
            path_cmap = clear_cmap(folder_out / f"{self.ms.molname}.{title}.cmap")
            vg.GridIO.write_cmap(path_cmap, self, self.ms.molname)
            return

        if gf == vg.GridFormat.CMAP_PACKED:
            path_cmap = clear_cmap(folder_out / f"{self.ms.molname}.cmap")
            vg.REMOVE_OLD_CMAP_OUTPUT = False
            vg.GridIO.write_cmap(path_cmap, self, f"{self.ms.molname}.{title}")
            return

        raise ValueError(f"Unknown output format: {vg.GRID_FORMAT_OUTPUT}.")


# //////////////////////////////////////////////////////////////////////////////
