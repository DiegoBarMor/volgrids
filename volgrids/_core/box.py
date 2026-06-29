import numpy as np

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Box:
    def __init__(self,
        origin: np.ndarray, resolution: np.ndarray, deltas: np.ndarray,
        do_init = True
    ):
        self.min_coords : np.ndarray[float] # minimum coordinates of the bounding box
        self.max_coords : np.ndarray[float] # maximum coordinates of the bounding box
        self.resolution : np.ndarray[int]   # number of grid points in each dimension
        self.deltas     : np.ndarray[float] # spacing between grid point in each dimension (in armstrong)
        self.cog        : np.ndarray[float] # center of geometry of the bounding box
        self.radius     : float             # (maximum) radius of the bounding box
        self.info       : vg.BoxInfo

        if not do_init: return # useful if values are to be set immediately after the Box initialization

        min_coords = np.array(origin, dtype = vg.FLOAT_DTYPE) # ensure they are numpy arrays...
        resolution = np.array(resolution, dtype = int)
        deltas     = np.array(deltas, dtype = vg.FLOAT_DTYPE)

        self.min_coords = min_coords
        self.max_coords = (min_coords + deltas * resolution).astype(vg.FLOAT_DTYPE)
        self.resolution = resolution
        self.deltas     = deltas
        self.infer_radius_and_cog()

        self.info = vg.BoxInfo.from_box(self)


    # --------------------------------------------------------------------------
    def __repr__(self) -> str:
        return f"Box({self.min_coords}->{self.max_coords})"


    # --------------------------------------------------------------------------
    @classmethod
    def from_min_max(cls, min_coords: np.ndarray, max_coords: np.ndarray) -> "Box":
        if np.any(max_coords <= min_coords):
            xmin, ymin, zmin = min_coords
            xmax, ymax, zmax = max_coords
            raise ValueError(
                f"Can't create box with min=({xmin},{ymin},{zmin}) and max=({xmax},{ymax},{zmax}): "
                f"Each box max coordinate must be strictly greater than its min."
            )

        box = cls(None, None, None, do_init = False)
        box.min_coords = np.array(min_coords)
        box.max_coords = np.array(max_coords)
        box.infer_deltas_resolution()
        box.infer_radius_and_cog()
        box.info = vg.BoxInfo.from_box(box)
        return box


    # --------------------------------------------------------------------------
    def __eq__(self, other: "Box") -> bool:
        return all((
            np.allclose(self.min_coords, other.min_coords),
            np.allclose(self.max_coords, other.max_coords),
            np.array_equal(self.resolution, other.resolution),
            np.allclose(self.deltas, other.deltas),
        ))


    # --------------------------------------------------------------------------
    def __neq__(self, other: "Box") -> bool:
        return not self.__eq__(other)


    # --------------------------------------------------------------------------
    def infer_deltas_resolution(self):
        box_size: np.ndarray = self.max_coords - self.min_coords
        if vg.CFG.box_from_deltas:
            self.deltas = np.array([vg.CFG.box_dx, vg.CFG.box_dy, vg.CFG.box_dz])
            self.resolution = np.round(box_size / self.deltas).astype(int)
        else:
            self.resolution = np.array([vg.CFG.box_xres, vg.CFG.box_yres, vg.CFG.box_zres], dtype = int)
            self.deltas = box_size / self.resolution


    # --------------------------------------------------------------------------
    def infer_radius_and_cog(self):
        self.radius = np.linalg.norm(self.max_coords - self.min_coords) / 2
        self.cog = (self.min_coords + self.max_coords) / 2


    # --------------------------------------------------------------------------
    def contains_point(self, point: np.ndarray) -> bool:
        return np.all((point >= self.min_coords) & (point <= self.max_coords))


    # --------------------------------------------------------------------------
    def pos_to_index(self, pos: np.ndarray[float]) -> np.ndarray[int]|None:
        """
        Convert a position in space to the corresponding index in the grid array.
        input:  (3,) FLOAT
        output: (3,) INT
            or `None` if the position is outside the box
        """
        if not self.contains_point(pos): return
        return np.round((pos - self.min_coords) / self.deltas).astype(int)


    # --------------------------------------------------------------------------
    def enforce_equilateral(self):
        max_resolution: np.ndarray = np.max(self.resolution)
        res_diff = max_resolution - self.resolution

        pad_0 = np.ceil (res_diff / 2).astype(int)
        pad_1 = np.floor(res_diff / 2).astype(int)

        self.resolution  = np.array([max_resolution, max_resolution, max_resolution], dtype = int)
        self.min_coords -= pad_0 * self.deltas
        self.max_coords += pad_1 * self.deltas
        self.infer_radius_and_cog()


    # --------------------------------------------------------------------------
    @classmethod
    def smallest_enclosing_box(cls, *boxes: "Box") -> "Box":
        return cls.from_min_max(
            min_coords = np.min([box.min_coords for box in boxes], axis = 0),
            max_coords = np.max([box.max_coords for box in boxes], axis = 0),
        )


# //////////////////////////////////////////////////////////////////////////////
