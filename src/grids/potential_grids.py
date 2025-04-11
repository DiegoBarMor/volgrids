import numpy as np

from src.utilities import Timer, write_json

from settings import DO_DX_OUTPUT, DO_MRC_OUTPUT, DO_JSON_OUTPUT, VERBOSE


# //////////////////////////////////////////////////////////////////////////////
class PotentialGrid:
    POTENTIAL_TYPE = "null" # OVERRIDE THIS

    def __init__(self, ms):
        if ms is None: return

        timer = Timer(f"...>>> PG: Creating {self.POTENTIAL_TYPE} potential grid...")
        self.ms = ms
        self.grid = np.zeros(self.ms.resolution, dtype = np.float32)
        self.populate_grid()
        self.grid[self.ms.trimming_mask] = 0 # apply precalculated trimming mask
        timer.end()

        if VERBOSE: timer = Timer(f"......       Saving {self.POTENTIAL_TYPE} potential grid...")
        self.pack_data()
        self.save_data()
        if VERBOSE: timer.end()


    @classmethod
    def grid_diff(cls, pg_0, pg_1, name = "diff"):
        pg = cls(None)
        pg.ms = pg_0.ms
        pg.grid = pg_0.grid - pg_1.grid

        timer = Timer(f".../// PG: Saving {name} potential grid...")
        pg.pack_data()
        pg.save_data(name)
        timer.end()


    def pack_data(self):
        xres, yres, zres = self.ms.resolution
        xmin, ymin, zmin = self.ms.minCoords
        xdelta, ydelta, zdelta = self.ms.deltas
        self.data = {
            "name" : self.ms.name,
            "xResolution" : int(xres - 1),
            "yResolution" : int(yres - 1),
            "zResolution" : int(zres - 1),
            "xOrigin" : float(xmin),
            "yOrigin" : float(ymin),
            "zOrigin" : float(zmin),
            "xDelta" : float(xdelta),
            "yDelta" : float(ydelta),
            "zDelta" : float(zdelta),
            "radius" : float(self.ms.radius),

            "minPotential" : float(np.min(self.grid)),
            "maxPotential" : float(np.max(self.grid)),
            "potentials" : self.grid.T.astype(float).flatten().tolist(),
        }


    def save_data(self, override_prefix = None):
        prefix = self.POTENTIAL_TYPE if override_prefix is None else override_prefix
        path_prefix = self.ms.folder_potentials / f"{self.ms.name}.{prefix}"

        ###### save to dx
        if DO_DX_OUTPUT: self.ms.write_dx(f"{path_prefix}.dx", self.grid)

        ###### save to mrc
        if DO_MRC_OUTPUT: self.ms.write_mrc(f"{path_prefix}.mrc", self.grid)

        ###### save to json
        if DO_JSON_OUTPUT: write_json(f"{path_prefix}.json", self.data)


    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
    def populate_grid(self): return


# //////////////////////////////////////////////////////////////////////////////
class StatisticalPotentialGrid(PotentialGrid):
    def populate_grid(self):
        for particle in self.iter_particles():
            self.process_particle(particle)

    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
    def iter_particles(self): return
    def process_particle(self, particle): return


# //////////////////////////////////////////////////////////////////////////////
