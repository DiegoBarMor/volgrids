import numpy as np
import volgrids as vg
from pathlib import Path
import MDAnalysis as mda

# //////////////////////////////////////////////////////////////////////////////
class MolecularSystem:
    def __init__(self, metadata):
        self.resolution = np.zeros(3, dtype = int)
        self.minCoords = np.zeros(3)
        self.deltas = np.zeros(3)

        self.path_pdb = Path(metadata["pdb"])
        self.path_apbs = Path(metadata["apbs"])
        self.path_metadata = Path(metadata["meta"])
        self.folder_potentials = Path(metadata["out"])

        self.name = self.path_pdb.stem
        self.isNucleic = bool(metadata["rna"])
        self.macro_query = f"{'nucleic' if self.isNucleic else 'protein'} and not (name H*)"

        if str(metadata["traj"]) == '.':
            self.system = mda.Universe(metadata["pdb"])
        else:
            self.system = mda.Universe(metadata["pdb"], metadata["traj"])

        self.frame = 0
        self.metadata = metadata.copy()

        ###############################
        self.radius : np.array
        self.cog : np.array
        self.minCoords : np.array
        self.maxCoords : np.array
        self.resolution : np.array
        self.deltas : np.array
        self.set_box_properties()


        ###############################
        grid_size = np.prod(self.resolution)
        if grid_size > vg.WARNING_GRID_SIZE:
            print()
            rx, ry, rz = self.resolution
            while True:
                choice = input(f">>> WARNING: resulting ({rx}x{ry}x{rz}) grid would contain {grid_size/1e6:.2f} million points. Proceed? [Y/N]\n").upper()
                if choice == 'Y': break
                if choice == 'N': exit()

    @classmethod
    def simple_init(cls, path_pdb: Path, folder_out: Path):
        return cls(metadata = {
            "pdb": path_pdb, "out": folder_out, "apbs": '',
            "rna": False, "default_res": False,
            "meta": folder_out / f"{path_pdb.stem}.meta.json",
        })

    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
    def set_box_properties(self): return
    def get_relevant_atoms(self): return
    def get_relevant_atoms_broad(self, trimming_dist): return


# //////////////////////////////////////////////////////////////////////////////
class MSPocketSphere(MolecularSystem):
    def set_box_properties(self):
        self.cog = np.array([self.metadata["xcog"], self.metadata["ycog"], self.metadata["zcog"]])
        self.radius = self.metadata["radius"]
        self.minCoords = self.cog - self.radius
        self.maxCoords = self.cog + self.radius
        box_size = self.maxCoords - self.minCoords

        if self.metadata["default_res"]:
            self.resolution = np.array([vg.GRID_XRES_PS, vg.GRID_YRES_PS, vg.GRID_ZRES_PS])
            self.deltas = box_size / self.resolution
        else:
            self.deltas = np.array([vg.GRID_DX_PS, vg.GRID_DY_PS, vg.GRID_DZ_PS])
            self.resolution = np.round(box_size / self.deltas).astype(int)


    def get_relevant_atoms(self):
        xcog, ycog, zcog = self.cog
        return self.system.select_atoms(
            f"{self.macro_query} and point {xcog} {ycog} {zcog} {self.radius}"
        )


    def get_relevant_atoms_broad(self, trimming_dist):
        xcog, ycog, zcog = self.cog
        return self.system.select_atoms(
            f"{self.macro_query} and point {xcog} {ycog} {zcog} {self.radius + trimming_dist}"
        )


# //////////////////////////////////////////////////////////////////////////////
class MSWhole(MolecularSystem):
    def set_box_properties(self):
        self.maxCoords = np.max(self.system.coord.positions, axis = 0) + vg.EXTRA_BOX_SIZE
        self.minCoords = np.min(self.system.coord.positions, axis = 0) - vg.EXTRA_BOX_SIZE
        self.radius = np.linalg.norm(self.maxCoords - self.minCoords) / 2
        self.cog = (self.maxCoords + self.minCoords) / 2
        box_size = self.maxCoords - self.minCoords

        if self.metadata["default_res"]:
            self.resolution = np.array([vg.GRID_XRES_WM, vg.GRID_YRES_WM, vg.GRID_ZRES_WM])
            self.deltas = box_size / self.resolution
        else:
            self.deltas = np.array([vg.GRID_DX_WM, vg.GRID_DY_WM, vg.GRID_DZ_WM])
            self.resolution = np.round(box_size / self.deltas).astype(int)


    def get_relevant_atoms(self):
        return self.system.select_atoms(self.macro_query)


    def get_relevant_atoms_broad(self, trimming_dist):
        return self.system.select_atoms(self.macro_query)


# //////////////////////////////////////////////////////////////////////////////
