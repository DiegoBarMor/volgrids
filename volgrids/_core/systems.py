import numpy as np
import MDAnalysis as mda
import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class MolecularSystem:
    def __init__(self, meta: "vg.ArgsParser"):
        self.resolution: np.array = np.zeros(3, dtype = int)
        self.minCoords : np.array = np.zeros(3)
        self.deltas    : np.array = np.zeros(3)
        self.meta: "vg.ArgsParser" = meta


        ###############################
        if meta.do_traj:
            self.system = mda.Universe(meta.path_structure, meta.path_traj)
            self.frame = 0
        else:
            self.system = mda.Universe(meta.path_structure)
            self.frame = None


        ###############################
        self.minCoords: np.array
        self.maxCoords: np.array
        self.radius: float
        self.cog: np.array
        self._init_box_attributes()


        ###############################
        box_size = self.maxCoords - self.minCoords
        if vg.USE_FIXED_DELTAS:
            self.deltas = np.array([vg.GRID_DX, vg.GRID_DY, vg.GRID_DZ])
            self.resolution = np.round(box_size / self.deltas).astype(int)
        else:
            self.resolution = np.array([vg.GRID_XRES, vg.GRID_YRES, vg.GRID_ZRES])
            self.deltas = box_size / self.resolution


        ###############################
        rx, ry, rz = self.resolution
        grid_size = rx*ry*rz
        if grid_size > vg.WARNING_GRID_SIZE:
            print()
            while True:
                choice = input(f">>> WARNING: resulting ({rx}x{ry}x{rz}) grid would contain {grid_size/1e6:.2f} million points. Proceed? [Y/N]\n").upper()
                if choice == 'Y': break
                if choice == 'N': exit()

    # --------------------------------------------------------------------------
    def _init_box_attributes(self):
        self.minCoords = np.min(self.system.coord.positions, axis = 0) - vg.EXTRA_BOX_SIZE
        self.maxCoords = np.max(self.system.coord.positions, axis = 0) + vg.EXTRA_BOX_SIZE
        self.radius = np.linalg.norm(self.maxCoords - self.minCoords) / 2
        self.cog = (self.minCoords + self.maxCoords) / 2

    # @classmethod
    # def simple_init(cls, path_pdb: Path, folder_out: Path):
    #     return cls(metadata = {
    #         "pdb": path_pdb, "out": folder_out, "apbs": '',
    #         "rna": False,
    #     })


# //////////////////////////////////////////////////////////////////////////////
