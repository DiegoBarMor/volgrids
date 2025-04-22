import numpy as np
import volpot as vp
from gridData import mrc

# //////////////////////////////////////////////////////////////////////////////
class Grid:
    def __init__(self, ms, dtype = np.float32):
        self.ms = ms
        self.xres, self.yres, self.zres = ms.resolution
        self.xmin, self.ymin, self.zmin = ms.minCoords
        self.dx, self.dy, self.dz = ms.deltas
        self.grid = np.zeros(ms.resolution, dtype = dtype)

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
        with mrc.mrcfile.new(path_mrc, overwrite = True) as mrc_file:
            mrc_file.set_data(self.grid.T.astype(np.float32))
            mrc_file.voxel_size = [self.dx, self.dy, self.dz]
            mrc_file.header["origin"]['x'] = self.xmin
            mrc_file.header["origin"]['y'] = self.ymin
            mrc_file.header["origin"]['z'] = self.zmin
            mrc_file.update_header_from_data()
            mrc_file.update_header_stats()

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
            f"object 3 class array type {dtype} rank 0 items {np.prod(self.ms.resolution)} data follows",
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
