import numpy as np
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Grid:
    def __init__(self, box: "vg.Box", init_grid = True, dtype = None):
        self.box = box
        self.dtype: type = vg.FLOAT_DTYPE if dtype is None else dtype
        self.arr: np.ndarray|None = np.zeros(box.resolution, dtype = self.dtype) if init_grid else None
        self.fmt: vg.GridFormat = None


    # --------------------------------------------------------------------------
    def __add__(self, other: "Grid|float|int") -> "Grid":
        obj = self.__class__(self.box, init_grid = False)
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
        obj = self.__class__(self.box, init_grid = False)
        if isinstance(other, Grid):
            obj.arr = self.arr - other.arr
            return obj
        try:
            obj.arr = self.arr - other
            return obj
        except TypeError:
            raise TypeError(f"Cannot substract {type(other)} from Grid. Use another Grid or a numeric value.")


    # --------------------------------------------------------------------------
    def __mul__(self, other: "Grid|float|int") -> "Grid":
        obj = Grid(self.ms, init_grid = False)
        if isinstance(other, Grid):
            obj.arr = self.arr * other.arr
            return obj
        try:
            obj.arr = self.arr * other
            return obj
        except TypeError:
            raise TypeError(f"Cannot multiply {type(other)} with Grid. Use another Grid or a numeric value.")


    # --------------------------------------------------------------------------
    def __abs__(self) -> "Grid":
        obj = self.__class__(self.box, init_grid = False)
        obj.arr = np.abs(self.arr)
        return obj


    # --------------------------------------------------------------------------
    @classmethod
    def reverse(cls, other: "Grid") -> "Grid":
        """Return a new Grid with the reversed values of the other Grid.
        For boolean grids, the reverse is the logical not.
        For numeric grids, the reverse is the negation of the values.
        """
        obj = cls(other.box, init_grid = False)
        obj.arr = np.logical_not(other.arr) if (other.dtype == bool) else -other.arr
        return obj


    # --------------------------------------------------------------------------
    def copy(self):
        obj = self.__class__(self.box, init_grid = False)
        obj.arr = np.copy(self.arr)
        return obj


    # -------------------------------------------------------------------------- GETTERS
    def xres(self): return self.box.resolution[0]
    def yres(self): return self.box.resolution[1]
    def zres(self): return self.box.resolution[2]
    def xmin(self): return self.box.min_coords[0]
    def ymin(self): return self.box.min_coords[1]
    def zmin(self): return self.box.min_coords[2]
    def xmax(self): return self.box.max_coords[0]
    def ymax(self): return self.box.max_coords[1]
    def zmax(self): return self.box.max_coords[2]
    def   dx(self): return self.box.deltas[0]
    def   dy(self): return self.box.deltas[1]
    def   dz(self): return self.box.deltas[2]

    def npoints(self): return self.xres() * self.yres() * self.zres()


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

        self.box.min_coords = np.array([new_xmin, new_ymin, new_zmin])
        self.box.max_coords = np.array([new_xmax, new_ymax, new_zmax])
        self.box.resolution = np.array([new_xres, new_yres, new_zres])
        self.box.deltas = (self.box.max_coords - self.box.min_coords) / (self.box.resolution - 1)
        self.box.infer_radius_and_cog()


    # --------------------------------------------------------------------------
    def has_equivalent_box(self, box: "vg.Box") -> bool:
        return not np.any(
            (self.box.resolution - box.resolution) +\
            (self.box.min_coords - box.min_coords) +\
            (self.box.max_coords - box.max_coords)
        )


    # --------------------------------------------------------------------------
    def reshape_as_box(self, box: "vg.Box"):
        self.reshape(
            new_min = box.min_coords,
            new_max = box.max_coords,
            new_res = box.resolution,
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
