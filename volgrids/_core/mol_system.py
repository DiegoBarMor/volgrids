import numpy as np
import MDAnalysis as mda
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class MolSystem:
    def __init__(self,
        path_struct: Path = None, path_traj: Path = None,
        box_data: dict = None
    ):
        self.minCoords  : np.ndarray[float]   # minimum coordinates of the bounding box
        self.maxCoords  : np.ndarray[float]   # maximum coordinates of the bounding box
        self.resolution : np.ndarray[int]     # number of grid points in each dimension
        self.deltas     : np.ndarray[float]   # size of each grid point in each dimension
        self.cog        : np.ndarray[float]   # center of geometry of the bounding box
        self.radius     : float               # (maximum) radius of the bounding box
        self.molname    : str                 # name of the molecule
        self.do_traj    : None | bool         # whether this is a trajectory or a single structure (None if no structure is provided)
        self.system     : None | mda.Universe # MDAnalysis Universe object for the molecular system
        self.frame      : None | int          # current frame number (if trajectory is used)

        if path_struct is not None:
            ### molecular system with a molecular structure (optionally a trajectory)
            ### the bounding box values are calculated from the structure
            self._init_attrs_from_molecules(path_struct, path_traj)

        else:
            ### simple molecular system with no molecular structure
            ### the bounding box values are to be expected instead
            if box_data is None:
                raise ValueError("Either 'path_struct' or 'box_data' must be provided.")
            self._init_attrs_from_box_data(box_data)


    # --------------------------------------------------------------------------
    @classmethod
    def from_box_data(cls,
        resolution: np.ndarray, origin: np.ndarray, deltas: np.ndarray,
        molname: str = "grid"
    ) -> "MolSystem":
        """
        Create a MolSystem instance from box data. 'maxCoords' is inferred
        :param resolution: The resolution of the grid (number of points in each dimension).
        :param origin: The origin of the grid (minimum coordinates of the bounding box).
        :param deltas: The size of each grid point in each dimension.
        :param molname: The name of the molecule (default is "grid").
        :return: An instance of MolSystem with the provided box data.
        """
        return cls(box_data = {
            "resolution": resolution,
            "minCoords": origin,
            "maxCoords": origin + deltas * resolution,
            "deltas": deltas,
            "molname": molname
        })


    # --------------------------------------------------------------------------
    def _init_attrs_from_molecules(self, path_struct: Path, path_traj: Path = None):
        self.molname = path_struct.stem
        self.do_traj = path_traj is not None

        if self.do_traj:
            self.system = mda.Universe(str(path_struct), str(path_traj))
            self.frame = 0
        else:
            self.system = mda.Universe(str(path_struct))
            self.frame = None

        self._infer_box_attributes()

        box_size = self.maxCoords - self.minCoords
        if vg.USE_FIXED_DELTAS:
            self.deltas = np.array([vg.GRID_DX, vg.GRID_DY, vg.GRID_DZ])
            self.resolution = np.round(box_size / self.deltas).astype(int)
        else:
            self.resolution = np.array([vg.GRID_XRES, vg.GRID_YRES, vg.GRID_ZRES], dtype = int)
            self.deltas = box_size / self.resolution

        self._warning_big_grid()


    # --------------------------------------------------------------------------
    def _init_attrs_from_box_data(self, box_data: dict):
        keys_box_data = set(box_data.keys())
        required_keys = {"molname", "minCoords", "maxCoords", "resolution", "deltas"}
        if not keys_box_data.issuperset(required_keys):
            raise ValueError(f"Box data must contain the keys: {required_keys}. Provided keys: {keys_box_data}")

        self.molname = box_data["molname"]
        self.do_traj = False

        self.system  = None
        self.frame   = None

        self.minCoords  = np.array(box_data["minCoords"],  dtype = float)
        self.maxCoords  = np.array(box_data["maxCoords"],  dtype = float)
        self.resolution = np.array(box_data["resolution"], dtype = int  )
        self.deltas     = np.array(box_data["deltas"],     dtype = float)
        self._calc_radius_and_cog()

        self._warning_big_grid()


    # --------------------------------------------------------------------------
    def _infer_box_attributes(self):
        self.minCoords = np.min(self.system.coord.positions, axis = 0) - vg.EXTRA_BOX_SIZE
        self.maxCoords = np.max(self.system.coord.positions, axis = 0) + vg.EXTRA_BOX_SIZE
        self._calc_radius_and_cog()


    # --------------------------------------------------------------------------
    def _calc_radius_and_cog(self):
        self.radius = np.linalg.norm(self.maxCoords - self.minCoords) / 2
        self.cog = (self.minCoords + self.maxCoords) / 2


    # --------------------------------------------------------------------------
    def _warning_big_grid(self):
        rx, ry, rz = self.resolution
        grid_size = rx*ry*rz
        if grid_size > vg.WARNING_GRID_SIZE:
            print()
            while True:
                choice = input(f">>> WARNING: resulting ({rx}x{ry}x{rz}) grid would contain {grid_size/1e6:.2f} million points. Proceed? [Y/N]\n").upper()
                if choice == 'Y': break
                if choice == 'N': exit()


# //////////////////////////////////////////////////////////////////////////////
