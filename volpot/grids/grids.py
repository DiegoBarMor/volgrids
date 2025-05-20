import numpy as np
import volpot as vp
import gridData as gd

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

    @classmethod
    def from_mrc(cls, path_mrc):
        with gd.mrc.mrcfile.open(path_mrc) as grid_mrc:
            vsize = np.array(grid_mrc.voxel_size)
            origin = np.array(grid_mrc.header["origin"])

            delta = np.array([vsize["x"], vsize["y"], vsize["z"]], dtype = np.float32)
            resolution = np.array(grid_mrc.data.shape)
            minCoords = np.array([origin["x"], origin["y"], origin["z"]], dtype = np.float32)

            grid_vp = cls(
                data = {
                    "resolution": resolution,
                    "minCoords": minCoords,
                    "maxCoords": minCoords + delta * resolution,
                    "deltas": delta
                },
                init_grid = False
            )
            grid_vp.grid = grid_mrc.data.T
        return grid_vp


    @classmethod
    def from_dx(cls, path_dx):
        grid_dx = gd.Grid(path_dx)
        grid_vp = cls(
            data = {
                "resolution": grid_dx.grid.shape,
                "minCoords": grid_dx.origin,
                "maxCoords": grid_dx.origin + grid_dx.delta * grid_dx.grid.shape,
                "deltas": grid_dx.delta
            },
            init_grid = False
        )
        grid_vp.grid = grid_dx.grid
        return grid_vp


    def run(self, trimmer: "vp.GridTrimmer"= None):
        timer = vp.Timer(f"...>>> PG: Calculating {self.get_type()} potential grid...")

        self.populate_grid()
        if trimmer is not None:
            trimmer.apply_trimming(self)

        self.data = {}
        self.pack_data()
        self.save_data()
        timer.end()

    def get_type(self):
        return ### OVERRIDE

    def copy(self):
        pg = Grid(self.ms)
        pg.grid = np.copy(self.grid)
        return pg

    def is_empty(self):
        return np.all(self.grid == 0)

    @classmethod
    def grid_diff(cls, pg_0, pg_1, name = "diff"):
        timer = vp.Timer(f".../// PG: Calculating {name} potential grid...")
        pg = pg_0.copy()
        pg.grid = pg_0.grid - pg_1.grid
        pg.pack_data()
        pg.save_data(name)
        timer.end()


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


    def save_data(self, override_prefix = None):
        prefix = self.get_type() if override_prefix is None else override_prefix
        path_prefix = self.ms.folder_potentials / f"{self.ms.name}.{prefix}"

        ###### save to dx
        if vp.DO_DX_OUTPUT: self.write_dx(f"{path_prefix}.dx")

        ###### save to mrc
        if vp.DO_MRC_OUTPUT: self.write_mrc(f"{path_prefix}.mrc")

        ###### save to json
        if vp.DO_JSON_OUTPUT: vp.write_json(f"{path_prefix}.json", self.data)


    # --------------------------------------------------------------------------
    def write_mrc(self, path_mrc):
        with gd.mrc.mrcfile.new(path_mrc, overwrite = True) as grid_mrc:
            grid_mrc.set_data(self.grid.T.astype(np.float32))
            grid_mrc.voxel_size = [self.dx, self.dy, self.dz]
            grid_mrc.header["origin"]['x'] = self.xmin
            grid_mrc.header["origin"]['y'] = self.ymin
            grid_mrc.header["origin"]['z'] = self.zmin
            grid_mrc.update_header_from_data()
            grid_mrc.update_header_stats()

    # ------------------------------------------------------------------------------
    def write_dx(self, path_dx):
        if self.grid.dtype == np.float32:
            grid_data = self.grid
            dtype = '"float"'
            fmt = "%.3f"
        elif self.grid.dtype == int:
            grid_data = self.grid
            dtype = '"int"'
            fmt = "%i"
        elif self.grid.dtype == bool:
            grid_data = self.grid.astype(int)
            dtype = '"int"'
            fmt = "%i"
        else:
            raise TypeError(f"Unsupported data type for DX output: {self.grid.dtype}")

        header = '\n'.join((
            "# OpenDX density file written by volpot2",
            "# File format: http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF",
            "# Data are embedded in the header and tied to the grid positions.",
            "# Data is written in C array order: In grid[x,y,z] the axis z is fastest",
            "# varying, then y, then finally x, i.e. z is the innermost loop.",
            f"object 1 class gridpositions counts {self.xres} {self.yres} {self.zres}",
            f"origin {self.xmin:6e} {self.ymin:6e} {self.zmin:6e}",
            f"delta {self.dx:6e} {0:6e} {0:6e}",
            f"delta {0:6e} {self.dy:6e} {0:6e}",
            f"delta {0:6e} {0:6e} {self.dz:6e}",
            f"object 2 class gridconnections counts  {self.xres} {self.yres} {self.zres}",
            f"object 3 class array type {dtype} rank 0 items {self.xres*self.yres*self.zres}, data follows",
        ))
        footer = '\n'.join((
            '',
            'attribute "dep" string "positions"',
            'object "density" class field',
            'component "positions" value 1',
            'component "connections" value 2',
            'component "data" value 3',
        ))

        ########### reshape the grid array
        grid_size = np.prod(grid_data.shape)
        dx_rows = grid_size // 3

        truncated_arr, extra_arr = np.split(grid_data.flatten(), [3*dx_rows])
        data_out = truncated_arr.reshape(dx_rows, 3)
        last_row = extra_arr.reshape(1, len(extra_arr))

        ########### export reshaped data
        with open(path_dx, "wb") as file:
            np.savetxt(
                file, data_out, fmt = fmt, delimiter = '\t',
                header = header, comments = ''
            )
            np.savetxt(
                file, last_row, fmt = fmt, delimiter = '\t',
                footer = footer, comments = ''
            )



    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
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
