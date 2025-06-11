import numpy as np
import MDAnalysis as mda
import volgrids.vgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class MolecularSystem:
    def __init__(self, meta: "vg.ArgsParser"):
        self.resolution: np.array = np.zeros(3, dtype = int)
        self.minCoords : np.array = np.zeros(3)
        self.deltas    : np.array = np.zeros(3)
        self.meta: "vg.ArgsParser" = meta

        self.isNucleic = meta.mode == "rna" # [WIP]
        self.macro_query = f"{meta.moltype} and not (name H*)" # [WIP]

        args =  (meta.path_in, meta.path_traj) if meta.do_traj else (meta.path_in, )
        self.system = mda.Universe(*args)

        self.frame = 0

        ###############################
        self.radius : np.array
        self.cog : np.array
        self.minCoords : np.array
        self.maxCoords : np.array
        self.set_box_properties()

        box_size = self.maxCoords - self.minCoords
        if vg.USE_FIXED_DELTAS:
            self.deltas = np.array([vg.GRID_DX, vg.GRID_DY, vg.GRID_DZ])
            self.resolution = np.round(box_size / self.deltas).astype(int)
        else:
            self.resolution = np.array([vg.GRID_XRES, vg.GRID_YRES, vg.GRID_ZRES])
            self.deltas = box_size / self.resolution


        ###############################
        grid_size = np.prod(self.resolution)
        if grid_size > vg.WARNING_GRID_SIZE:
            print()
            rx, ry, rz = self.resolution
            while True:
                choice = input(f">>> WARNING: resulting ({rx}x{ry}x{rz}) grid would contain {grid_size/1e6:.2f} million points. Proceed? [Y/N]\n").upper()
                if choice == 'Y': break
                if choice == 'N': exit()

    # @classmethod
    # def simple_init(cls, path_pdb: Path, folder_out: Path):
    #     return cls(metadata = {
    #         "pdb": path_pdb, "out": folder_out, "apbs": '',
    #         "rna": False,
    #         "meta": folder_out / f"{path_pdb.stem}.meta.json",
    #     })

    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
    def set_box_properties(self): return
    def get_relevant_atoms(self): return
    def get_relevant_atoms_broad(self, trimming_dist): return


# //////////////////////////////////////////////////////////////////////////////
class MSPocketSphere(MolecularSystem):
    def set_box_properties(self):
        r,x,y,z = self.meta.ps_info
        self.cog = np.array([x, y, z])
        self.radius = r
        self.minCoords = self.cog - self.radius
        self.maxCoords = self.cog + self.radius


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


    def get_relevant_atoms(self):
        return self.system.select_atoms(self.macro_query)


    def get_relevant_atoms_broad(self, trimming_dist):
        return self.system.select_atoms(self.macro_query)


# //////////////////////////////////////////////////////////////////////////////
