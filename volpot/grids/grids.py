import numpy as np
import volpot as vp

# //////////////////////////////////////////////////////////////////////////////
class Grid:
    def __init__(self, data, dtype = np.float32, init_grid = True):
        if isinstance(data, vp.MolecularSystem):
            self.ms = data
            self.xres, self.yres, self.zres = data.resolution
            self.xmin, self.ymin, self.zmin = data.minCoords
            self.xmax, self.ymax, self.zmax = data.maxCoords
            self.dx, self.dy, self.dz = data.deltas

        elif isinstance(data, dict):
            self.ms = None
            self.xres, self.yres, self.zres = data["resolution"]
            self.xmin, self.ymin, self.zmin = data["minCoords"]
            self.xmax, self.ymax, self.zmax = data["maxCoords"]
            self.dx, self.dy, self.dz = data["deltas"]

        self.grid = np.zeros((self.xres, self.yres, self.zres), dtype = dtype) \
            if init_grid else None

    # --------------------------------------------------------------------------
    def run(self, trimmer: "vp.GridTrimmer"= None):
        self.populate_grid()
        if trimmer is not None:
            trimmer.apply_trimming(self)

        self.data = {}
        self.pack_data()
        self.save_data()

    # --------------------------------------------------------------------------
    def copy(self):
        pg = Grid(self.ms)
        pg.grid = np.copy(self.grid)
        return pg

    # --------------------------------------------------------------------------
    def is_empty(self):
        return np.all(self.grid == 0)

    # --------------------------------------------------------------------------
    @classmethod
    def grid_diff(cls, pg_0, pg_1, name = "diff"):
        pg = pg_0.copy()
        pg.grid = pg_0.grid - pg_1.grid
        pg.pack_data()
        pg.save_data(name)

    # --------------------------------------------------------------------------
    def pack_data(self):
        if vp.DO_JSON_OUTPUT:
            self.data = {
                "name" : self.ms.name,
                "xResolution" : int(self.xres - 1),
                "yResolution" : int(self.yres - 1),
                "zResolution" : int(self.zres - 1),
                "xOrigin" : float(self.xmin),
                "yOrigin" : float(self.ymin),
                "zOrigin" : float(self.zmin),
                "xDelta" : float(self.dx),
                "yDelta" : float(self.dy),
                "zDelta" : float(self.dz),
                "radius" : float(self.ms.radius),
                "minPotential" : float(np.min(self.grid)),
                "maxPotential" : float(np.max(self.grid)),
                "potentials" : self.grid.T.astype(float).flatten().tolist(), # [TODO] check if this should be transposed or not
            }

    # --------------------------------------------------------------------------
    def save_data(self, override_prefix = None):
        prefix = self.get_type() if override_prefix is None else override_prefix
        path_prefix = self.ms.folder_potentials / f"{self.ms.name}.{prefix}"

        ###### save to dx
        if vp.DO_DX_OUTPUT: vp.write_dx(f"{path_prefix}.dx", self)

        ###### save to mrc
        if vp.DO_MRC_OUTPUT: vp.write_mrc(f"{path_prefix}.mrc", self)

        ###### save to cmap
        if vp.DO_CMAP_OUTPUT: vp.write_cmap(f"{path_prefix}.cmap", self, f"{self.ms.name}.{self.ms.frame:04}")

        ###### save to json
        if vp.DO_JSON_OUTPUT: vp.write_json(f"{path_prefix}.json", self.data)


    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
    def get_type(self): return
    def populate_grid(self): return


# //////////////////////////////////////////////////////////////////////////////
class GridSMIF(Grid):
    def populate_grid(self):
        for particle in self.iter_particles():
            self.process_particle(particle)

    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
    def iter_particles(self): return
    def process_particle(self, particle): return


# //////////////////////////////////////////////////////////////////////////////
