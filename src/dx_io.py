import numpy as np
from gridData import mrc


# //////////////////////////////////////////////////////////////////////////////
class DXIO:
    def __init__(self):
        self.resolution = np.zeros(3, dtype = int)
        self.minCoords = np.zeros(3)
        self.deltas = np.zeros(3)

        self.dx_comments = '\n'.join((
            "# OpenDX density file written by dx_io.DXIO.write_dx()",
            "# File format: http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF",
            "# Data are embedded in the header and tied to the grid positions.",
            "# Data is written in C array order: In grid[x,y,z] the axis z is fastest",
            "# varying, then y, then finally x, i.e. z is the innermost loop.",
            '',
        ))
        self.dx_footer = '\n'.join((
            '',
            'attribute "dep" string "positions"',
            'object "density" class field',
            'component "positions" value 1',
            'component "connections" value 2',
            'component "data" value 3',
        ))


    # --------------------------------------------------------------------------
    def write_dx(self, path_dx, grid_data):
        ########### init values
        xres,yres,zres = self.resolution
        xmin,ymin,zmin = self.minCoords
        dx,dy,dz = self.deltas

        if grid_data.dtype == np.float32:
            dtype = '"float"'
            fmt = "%.3f"
        elif grid_data.dtype == int:
            dtype = '"int"'
            fmt = "%i"

        header = self.dx_comments + '\n'.join((
            f"object 1 class gridpositions counts  {xres} {yres} {zres}",
            f"origin {xmin:6f} {ymin:6f} {zmin:6f}",
            f"delta  {dx} 0 0",
            f"delta  0 {dy} 0",
            f"delta  0 0 {dz}",
            f"object 2 class gridconnections counts  {xres} {yres} {zres}",
            f"object 3 class array type {dtype} rank 0 items {np.prod(self.resolution)} data follows",
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
                footer = self.dx_footer, comments = ''
            )


    # --------------------------------------------------------------------------
    def write_mrc(self, path_mrc, grid_data):
        with mrc.mrcfile.new(path_mrc, overwrite = True) as mrc_file:
            mrc_file.set_data(grid_data.T.astype(np.float32))
            mrc_file.voxel_size = self.deltas.tolist()
            mrc_file.header["origin"]['x'] = self.minCoords[0]
            mrc_file.header["origin"]['y'] = self.minCoords[1]
            mrc_file.header["origin"]['z'] = self.minCoords[2]
            mrc_file.update_header_from_data()
            mrc_file.update_header_stats()


# //////////////////////////////////////////////////////////////////////////////
