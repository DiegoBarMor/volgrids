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

        if not do_init: return # useful if values are to be set immediately after the Box initialization

        min_coords = np.array(origin, dtype = vg.FLOAT_DTYPE) # ensure they are numpy arrays...
        resolution = np.array(resolution, dtype = int)
        deltas     = np.array(deltas, dtype = vg.FLOAT_DTYPE)

        self.min_coords = min_coords
        self.max_coords = (min_coords + deltas * resolution).astype(vg.FLOAT_DTYPE)
        self.resolution = resolution
        self.deltas     = deltas
        self.infer_radius_and_cog()

        self._warning_big_grid()


    # --------------------------------------------------------------------------
    def infer_deltas_resolution(self):
        box_size: np.ndarray = self.max_coords - self.min_coords
        if vg.USE_FIXED_DELTAS:
            self.deltas = np.array([vg.GRID_DX, vg.GRID_DY, vg.GRID_DZ])
            self.resolution = np.round(box_size / self.deltas).astype(int)
        else:
            self.resolution = np.array([vg.GRID_XRES, vg.GRID_YRES, vg.GRID_ZRES], dtype = int)
            self.deltas = box_size / self.resolution


    # --------------------------------------------------------------------------
    def infer_radius_and_cog(self):
        self.radius = np.linalg.norm(self.max_coords - self.min_coords) / 2
        self.cog = (self.min_coords + self.max_coords) / 2


    # --------------------------------------------------------------------------
    def _warning_big_grid(self):
        rx, ry, rz = self.resolution
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
