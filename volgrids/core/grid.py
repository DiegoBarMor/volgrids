import numpy as np
import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Grid:
    def __init__(self, data, init_grid = True, dtype = None):
        if isinstance(data, vg.MolecularSystem):
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

        if dtype is None: dtype = vg.FLOAT_DTYPE
        self.grid = np.zeros((self.xres, self.yres, self.zres), dtype = dtype) \
            if init_grid else None

    # --------------------------------------------------------------------------
    def run(self, trimmer: "vg.GridTrimmer"= None):
        self.populate_grid()
        if trimmer is not None:
            trimmer.apply_trimming(self)

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
    def reshape(self, new_xres, new_yres, new_zres):
        self.grid = vg.interpolate_3d(
            x0 = np.linspace(self.xmin, self.xmax, self.xres),
            y0 = np.linspace(self.ymin, self.ymax, self.yres),
            z0 = np.linspace(self.zmin, self.zmax, self.zres),
            data_0 = self.grid,
            new_coords = np.mgrid[
                self.xmin : self.xmax : complex(0, new_xres),
                self.ymin : self.ymax : complex(0, new_yres),
                self.zmin : self.zmax : complex(0, new_zres),
            ].T
        ).astype(vg.FLOAT_DTYPE)
        self.xres = new_xres
        self.yres = new_yres
        self.zres = new_zres

    # --------------------------------------------------------------------------
    @classmethod
    def grid_diff(cls, pg_0, pg_1, name = "diff"):
        pg = pg_0.copy()
        pg.grid = pg_0.grid - pg_1.grid
        pg.save_data(name)

    # --------------------------------------------------------------------------
    def save_data(self, override_prefix = None):
        name = self.ms.meta.name
        prefix = self.get_type() if override_prefix is None else override_prefix
        path_prefix = self.ms.meta.path_out / f"{name}.{prefix}"

        if self.ms.meta.do_traj:
            ### ignore the OUTPUT flags, CMAP is the only format that supports multiple frames
            vg.write_cmap(f"{path_prefix}.cmap", self, f"{name}.{self.ms.frame:04}")
            return

        if vg.DO_OUTPUT_DX:   vg.write_dx  (f"{path_prefix}.dx",   self)
        if vg.DO_OUTPUT_MRC:  vg.write_mrc (f"{path_prefix}.mrc",  self)
        if vg.DO_OUTPUT_CMAP: vg.write_cmap(f"{path_prefix}.cmap", self, name)


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
