import numpy as np
import volpot as vp
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

        self.system = mda.Universe(metadata["pdb"])
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
        if grid_size > vp.WARNING_GRID_SIZE:
            print()
            while True:
                choice = input(f">>> WARNING: resulting grid would contain {grid_size/1e6:.2f} million points. Proceed? [Y/N]\n").upper()
                if choice == 'Y': break
                if choice == 'N': exit()


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
            self.resolution = np.array([vp.GRID_XRES_PS, vp.GRID_YRES_PS, vp.GRID_ZRES_PS])
            self.deltas = box_size / self.resolution
        else:
            self.deltas = np.array([vp.GRID_DX_PS, vp.GRID_DY_PS, vp.GRID_DZ_PS])
            self.resolution = np.round(box_size / self.deltas).astype(int)

        print("...>>> MS: PocketSphere mode. Info:", flush = True)


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
        self.maxCoords = np.max(self.system.coord.positions, axis = 0) + vp.EXTRA_BOX_SIZE
        self.minCoords = np.min(self.system.coord.positions, axis = 0) - vp.EXTRA_BOX_SIZE
        self.radius = vp.get_norm(self.maxCoords - self.minCoords) / 2
        self.cog = (self.maxCoords + self.minCoords) / 2
        box_size = self.maxCoords - self.minCoords

        if self.metadata["default_res"]:
            self.resolution = np.array([vp.GRID_XRES_WM, vp.GRID_YRES_WM, vp.GRID_ZRES_WM])
            self.deltas = box_size / self.resolution
        else:
            self.deltas = np.array([vp.GRID_DX_WM, vp.GRID_DY_WM, vp.GRID_DZ_WM])
            self.resolution = np.round(box_size / self.deltas).astype(int)

        print("...>>> MS: Whole mode. Info:", flush = True)

    def get_relevant_atoms(self):
        return self.system.select_atoms(self.macro_query)

    def get_relevant_atoms_broad(self, trimming_dist):
        return self.system.select_atoms(self.macro_query)


# //////////////////////////////////////////////////////////////////////////////
