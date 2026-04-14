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
        obj = self.__class__(self.ms, init_grid = False)
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
        obj = self.__class__(self.ms, init_grid = False)
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
        obj = self.__class__(self.ms, init_grid = False)
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
        obj = self.__class__(self.ms, init_grid = False)
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
    def save_data(self, path_out: Path, grid_format: "vg.GridFormat", cmap_key = "grid"):
        path_out = Path(path_out)

        if grid_format == vg.GridFormat.DX:
            vg.GridIO.write_dx(path_out, self)
            return

        if grid_format == vg.GridFormat.MRC:
            vg.GridIO.write_mrc(path_out, self)
            return

        if grid_format == vg.GridFormat.CCP4:
            vg.GridIO.write_ccp4(path_out, self)
            return

        if grid_format == vg.GridFormat.CMAP:
            vg.GridIO.clear_cmap(path_out)
            vg.GridIO.write_cmap(path_out, self, cmap_key)
            return

        if grid_format == vg.GridFormat.CMAP_PACKED:
            vg.GridIO.clear_cmap(path_out)
            vg.REMOVE_OLD_CMAP_OUTPUT = False
            vg.GridIO.write_cmap(path_out, self, cmap_key)
            return

        raise ValueError(f"Unknown output format: {grid_format}.")


# //////////////////////////////////////////////////////////////////////////////
