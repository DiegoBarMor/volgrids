import numpy as np
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Grid:
    def __init__(self, box: "vg.Box", init_grid = True, dtype = None):
        self.box = box
        self.dtype: type = vg.FLOAT_DTYPE if dtype is None else dtype
        self.fmt: vg.GridFormat = None
        self.arr: np.ndarray|None
        self.dirty: bool = False # whether grid should be filled with 0s when calling reset()
        self.name: str = "" # mainly used as key when saving to cmap format

        if init_grid:
            self._warning_big_grid()
            self.arr = np.zeros(box.resolution, dtype = self.dtype)
        else:
            self.arr = None


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
        obj = self.__class__(self.box, init_grid = False)
        if isinstance(other, Grid):
            obj.arr = self.arr * other.arr
            return obj
        try:
            obj.arr = self.arr * other
            return obj
        except TypeError:
            raise TypeError(f"Cannot multiply {type(other)} with Grid. Use another Grid or a numeric value.")


    # --------------------------------------------------------------------------
    def __div__(self, other: "Grid|float|int") -> "Grid":
        obj = self.__class__(self.box, init_grid = False)
        if isinstance(other, Grid):
            obj.arr = self.arr / other.arr
            return obj
        try:
            obj.arr = self.arr / other
            return obj
        except TypeError:
            raise TypeError(f"Cannot divide Grid by {type(other)}. Use another Grid or a numeric value.")


    # --------------------------------------------------------------------------
    def __abs__(self) -> "Grid":
        obj = self.__class__(self.box, init_grid = False)
        obj.arr = np.abs(self.arr)
        return obj


    # --------------------------------------------------------------------------
    def __and__(self, other: "Grid") -> "Grid":
        """
        Return a new Grid with the logical or of the values of the two Grids.
        The input grids are converted into boolean if they are not already, and the output grid is also boolean.
        """
        obj = self.__class__(self.box, init_grid = False)
        if not isinstance(other, Grid):
            raise TypeError(f"Cannot add {type(other)} to Grid. Use another Grid.")
        vg.GridIO.restore_boolean_dtype(self)
        vg.GridIO.restore_boolean_dtype(other)
        obj.arr = np.logical_and(self.arr, other.arr)
        return obj


    # --------------------------------------------------------------------------
    def __or__(self, other: "Grid") -> "Grid":
        """
        Return a new Grid with the logical or of the values of the two Grids.
        The input grids are converted into boolean if they are not already, and the output grid is also boolean.
        """
        obj = self.__class__(self.box, init_grid = False)
        if not isinstance(other, Grid):
            raise TypeError(f"Cannot add {type(other)} to Grid. Use another Grid.")
        vg.GridIO.restore_boolean_dtype(self)
        vg.GridIO.restore_boolean_dtype(other)
        obj.arr = np.logical_or(self.arr, other.arr)
        return obj


    # --------------------------------------------------------------------------
    @classmethod
    def load(cls, path_in: Path, fmt: "vg.GridFormat" = None, key: str = None) -> "Grid":
        """
        If the format is not specified, it will be detected based on the input path extension.
        If the file is a CMAP file and key is not specified, it will read the first key found in the file.
        """
        path_in = Path(path_in)
        if fmt is None: fmt = vg.GridIO.detect_format(path_in)

        if fmt == vg.GridFormat.DX:
            return vg.GridIO.read_dx(path_in)

        if fmt == vg.GridFormat.BIN:
            return vg.GridIO.read_bin(path_in)

        if fmt == vg.GridFormat.MRC:
            return vg.GridIO.read_mrc(path_in)

        if fmt == vg.GridFormat.CCP4:
            return vg.GridIO.read_ccp4(path_in)

        if fmt == vg.GridFormat.CMAP:
            if key is None:
                keys = vg.GridIO.get_cmap_keys(path_in, assert_has_keys = True)
                key = keys[0]
            return vg.GridIO.read_cmap(path_in, key)

        raise ValueError(f"Unknown input format: {fmt}.")


    # --------------------------------------------------------------------------
    def save(self, path_out: Path, fmt: "vg.GridFormat" = None, key: str = ""):
        """
        If the format is not specified, it will be detected based on the output path extension.
        When saving to CMAP, the key passed via `key` will be prioritized over the `Grid.name` attribute.
        If neither is provided, a new key will be generated based on the number of existing keys in the CMAP file.
        """
        def _key_cmap():
            if key: return key
            if self.name: return self.name
            if not path_out.is_file(): return "grid_0000"
            keys = vg.GridIO.get_cmap_keys(path_out)
            return f"grid_{len(keys)+1:04d}"

        path_out = Path(path_out)
        if fmt is None: fmt = vg.GridIO.detect_format(path_out)

        if fmt == vg.GridFormat.DX:
            vg.GridIO.write_dx(path_out, self)
            return

        if fmt == vg.GridFormat.BIN:
            vg.GridIO.write_bin(path_out, self)
            return

        if fmt == vg.GridFormat.MRC:
            vg.GridIO.write_mrc(path_out, self)
            return

        if fmt == vg.GridFormat.CCP4:
            vg.GridIO.write_ccp4(path_out, self)
            return

        if fmt == vg.GridFormat.CMAP:
            vg.GridIO.write_cmap(path_out, self, _key_cmap())
            return

        raise ValueError(f"Unknown output format: {fmt}.")


    # --------------------------------------------------------------------------
    def reset(self, force = False) -> None:
        if force or self.dirty: self.arr.fill(0)


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
    def value_at_position(self, pos: tuple[float, float, float]) -> float:
        """Get the value of the grid at a given position in space. Returns `0.0` if the position is outside the grid."""
        idx = self.box.pos_to_index(np.array(pos))
        if idx is None: return 0.0
        return self.arr[tuple(idx)]


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
    def _warning_big_grid(self):
        rx, ry, rz = self.box.resolution
        grid_size = rx*ry*rz
        if grid_size < vg.WARNING_GRID_SIZE: return
        print()
        while True:
            choice = input(
                f">>> WARNING: resulting ({rx}x{ry}x{rz}) grid would contain {grid_size/1e6:.2f} million points. Proceed? [Y/N]\n"
            ).upper()
            if choice.startswith('Y'): break
            if choice.startswith('N'): exit(3)


# //////////////////////////////////////////////////////////////////////////////
