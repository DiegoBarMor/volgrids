import os, h5py
import numpy as np
import gridData as gd
from pathlib import Path
from enum import Enum, auto

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class GridFormat(Enum):
    DX = auto()
    MRC = auto()
    CCP4 = auto()
    CMAP = auto()
    CMAP_PACKED = auto()


# //////////////////////////////////////////////////////////////////////////////
class GridIO:
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN I/O OPERATIONS
    @staticmethod
    def read_dx(path_dx) -> "vg.Grid":
        parser_dx = gd.Grid(path_dx)
        ms = vg.MolSystem.from_box_data(
            resolution = parser_dx.grid.shape,
            origin = parser_dx.origin,
            deltas = parser_dx.delta
        )
        obj = vg.Grid(ms, init_grid = False)
        obj.grid = parser_dx.grid
        obj.fmt = vg.GridFormat.DX
        return obj


    # --------------------------------------------------------------------------
    @staticmethod
    def read_mrc(path_mrc) -> "vg.Grid":
        with gd.mrc.mrcfile.open(path_mrc) as parser:
            ##### assume that MRC always follows the origin follows the "real space" MRC convention
            orig = parser.header["origin"]
            used_origin = np.array([orig['x'], orig['y'], orig['z']])

        obj = _read_mrc_ccp4(path_mrc, used_origin)
        obj.fmt = vg.GridFormat.MRC
        return obj


    # --------------------------------------------------------------------------
    @staticmethod
    def read_ccp4(path_ccp4) -> "vg.Grid":
        with gd.mrc.mrcfile.open(path_ccp4) as parser:
            orig = parser.header["origin"]
            if (orig['x'] == 0.0 and orig['y'] == 0.0 and orig['z'] == 0.0):
                ##### assume the origin follows the "integer offset" CCP4 convention, so use that one
                ##### https://mail.cgl.ucsf.edu/mailman/archives/list/chimera-users@cgl.ucsf.edu/thread/GLP62Q2WIBJ6ZU4H6BWXWVXETYSXVIWS/
                used_origin = np.array([
                    parser.voxel_size['x'] * parser.header["nxstart"],
                    parser.voxel_size['y'] * parser.header["nystart"],
                    parser.voxel_size['z'] * parser.header["nzstart"],
                ])
            else:
                ##### assume the origin follows the "real space" MRC convention, so use that one
                used_origin = np.array([orig['x'], orig['y'], orig['z']])

        obj = _read_mrc_ccp4(path_ccp4, used_origin)
        obj.fmt = vg.GridFormat.CCP4
        return obj


    # --------------------------------------------------------------------------
    @staticmethod
    def read_cmap(path_cmap, key) -> "vg.Grid":
        with h5py.File(path_cmap, 'r') as parser:
            frame = parser["Chimera"][key]
            rz, ry, rx = frame["data_zyx"].shape
            ox, oy, oz = frame.attrs["origin"]
            dz, dy, dx = frame.attrs["step"]
            ms = vg.MolSystem.from_box_data(
                resolution = np.array([rx, ry, rz]),
                origin = np.array([ox, oy, oz]),
                deltas = np.array([dx, dy, dz])
            )
            obj = vg.Grid(ms, init_grid = False)
            obj.grid = frame["data_zyx"][()].transpose(2,1,0)

            n_keys = len(parser["Chimera"].keys())
            obj.fmt = vg.GridFormat.CMAP_PACKED \
                if (n_keys > 1) else vg.GridFormat.CMAP
        return obj


    # --------------------------------------------------------------------------
    @staticmethod
    def write_dx(path_dx, data: "vg.Grid"):
        ints = (int, np.int8, np.int16, np.int32, np.int64)
        floats = (float, np.float16, np.float32, np.float64)

        if data.grid.dtype in floats:
            grid_data = data.grid
            dtype = '"float"'
            fmt = "%.3f"
        elif data.grid.dtype in ints:
            grid_data = data.grid
            dtype = '"int"'
            fmt = "%i"
        elif data.grid.dtype == bool:
            grid_data = data.grid.astype(int)
            dtype = '"int"'
            fmt = "%i"
        else:
            raise TypeError(f"Unsupported data type for DX output: {data.grid.dtype}")

        header = '\n'.join((
            "# OpenDX density file written by volgrids",
            "# File format: http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF",
            "# Data are embedded in the header and tied to the grid positions.",
            "# Data is written in C array order: In grid[x,y,z] the axis z is fastest",
            "# varying, then y, then finally x, i.e. z is the innermost loop.",
            f"object 1 class gridpositions counts {data.xres} {data.yres} {data.zres}",
            f"origin {data.xmin:6e} {data.ymin:6e} {data.zmin:6e}",
            f"delta {data.dx:6e} {0:6e} {0:6e}",
            f"delta {0:6e} {data.dy:6e} {0:6e}",
            f"delta {0:6e} {0:6e} {data.dz:6e}",
            f"object 2 class gridconnections counts  {data.xres} {data.yres} {data.zres}",
            f"object 3 class array type {dtype} rank 0 items {data.xres*data.yres*data.zres}, data follows",
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


    # --------------------------------------------------------------------------
    @staticmethod
    def write_mrc(path_mrc, data: "vg.Grid"):
        with gd.mrc.mrcfile.new(path_mrc, overwrite = True) as parser:
            parser.set_data(data.grid.astype(vg.FLOAT_DTYPE).transpose(2,1,0))
            parser.voxel_size = [data.dx, data.dy, data.dz]
            parser.header["origin"]['x'] = data.xmin # MRC convention
            parser.header["origin"]['y'] = data.ymin
            parser.header["origin"]['z'] = data.zmin
            parser.update_header_from_data()
            parser.update_header_stats()


    # --------------------------------------------------------------------------
    @staticmethod
    def write_ccp4(path_ccp4, data: "vg.Grid"):
        with gd.mrc.mrcfile.new(path_ccp4, overwrite = True) as parser:
            parser.set_data(data.grid.astype(vg.FLOAT_DTYPE).transpose(2,1,0))
            parser.voxel_size = [data.dx, data.dy, data.dz]
            parser.header["origin"]['x'] = data.xmin # MRC convention
            parser.header["origin"]['y'] = data.ymin
            parser.header["origin"]['z'] = data.zmin
            parser.header["nxstart"] = int(data.xmin / data.dx) # CCP4 convention
            parser.header["nystart"] = int(data.ymin / data.dy)
            parser.header["nzstart"] = int(data.zmin / data.dz)
            parser.update_header_from_data()
            parser.update_header_stats()


    # --------------------------------------------------------------------------
    @staticmethod
    def write_cmap(path_cmap, data: "vg.Grid", key):
        ### imitate the Chimera cmap format, as "specified" in this sample:
        ### https://github.com/RBVI/ChimeraX/blob/develop/testdata/cell15_timeseries.cmap
        def _add_generic_attrs(group, c = "GROUP"):
            group.attrs["CLASS"] = np.bytes_(c)
            group.attrs["TITLE"] = np.bytes_("")
            group.attrs["VERSION"] = np.bytes_("1.0")

        if not os.path.exists(path_cmap):
            with h5py.File(path_cmap, 'w') as h5:
                h5.attrs["PYTABLES_FORMAT_VERSION"] = np.bytes_("2.0")
                _add_generic_attrs(h5)

                chim = h5.create_group("Chimera")
                _add_generic_attrs(chim)

        with h5py.File(path_cmap, 'a') as parser:
            chim = parser["Chimera"]
            if key in chim.keys():
                frame = chim[key]
                if "data_zyx" in frame.keys():
                    del frame["data_zyx"]
            else:
                frame = parser.create_group(f"/Chimera/{key}")
                frame.attrs["chimera_map_version"] = np.int64(1)
                frame.attrs["chimera_version"] = np.bytes_(b'1.12_b40875')
                frame.attrs["name"] = np.bytes_(key)
                frame.attrs["origin"] = np.array([data.xmin, data.ymin, data.zmin], dtype = vg.FLOAT_DTYPE)
                frame.attrs["step"] = np.array([data.dz, data.dy, data.dx], dtype = vg.FLOAT_DTYPE)
                _add_generic_attrs(frame)

            framedata = frame.create_dataset(
                "data_zyx", data = data.grid.transpose(2,1,0), dtype = vg.FLOAT_DTYPE,
                compression = "gzip", compression_opts = vg.GZIP_COMPRESSION
            )
            _add_generic_attrs(framedata, "CARRAY")


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ OTHER I/O UTILITIES
    @staticmethod
    def read_auto(path_grid: Path) -> "vg.Grid":
        """Detect the format of the grid file based on its extension and then read it."""
        ext = path_grid.suffix.lower()

        # [TODO] improve the format detection?
        if ext == ".dx":
            return GridIO.read_dx(path_grid)

        if ext == ".mrc":
            return GridIO.read_mrc(path_grid)

        if ext == ".ccp4":
            return GridIO.read_ccp4(path_grid)

        if ext == ".cmap":
            keys = GridIO.get_cmap_keys(path_grid)
            if not keys: raise ValueError(f"Empty cmap file: {path_grid}")
            return GridIO.read_cmap(path_grid, keys[0])

        raise ValueError(f"Unrecognized file format: {ext}")


    # --------------------------------------------------------------------------
    @staticmethod
    def write_auto(path_grid: Path, data: "vg.Grid"):
        """Detect the format of the grid file based on its extension and then write it."""
        ext = path_grid.suffix.lower()

        # [TODO] improve the format detection?
        if ext == ".dx":
            GridIO.write_dx(path_grid, data)
            return

        if ext == ".mrc":
            GridIO.write_mrc(path_grid, data)
            return

        if ext == ".ccp4":
            GridIO.write_ccp4(path_grid, data)
            return

        if ext == ".cmap":
            keys = GridIO.get_cmap_keys(path_grid)
            GridIO.write_cmap(path_grid, data, key = f"map_{len(keys)+1}")
            return

        raise ValueError(f"Unrecognized file format: {ext}")


    # --------------------------------------------------------------------------
    @staticmethod
    def get_cmap_keys(path_cmap) -> list[str]:
        with h5py.File(path_cmap, 'r') as h5:
            return list(h5["Chimera"].keys())


# //////////////////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------------------
def _read_mrc_ccp4(path_mrc, origin: np.ndarray) -> "vg.Grid":
    with gd.mrc.mrcfile.open(path_mrc) as parser:
        # machine_stamp = parser.header.machst
        ### [68 68 0 0] or [68 65 0 0] for little-endian <--- tested
        ### [17 17 0 0] for big-endian <--- what happens in these cases?

        vsize = np.array([
            parser.voxel_size['x'],
            parser.voxel_size['y'],
            parser.voxel_size['z'],
        ], dtype = vg.FLOAT_DTYPE)

        res = np.array([
            parser.header["mx"],
            parser.header["my"],
            parser.header["mz"],
        ], dtype = int)

        data: np.ndarray = parser.data.astype(vg.FLOAT_DTYPE)
        origin = origin.astype(vg.FLOAT_DTYPE)

        axes_correspondance =\
            parser.header.mapc, parser.header.mapr, parser.header.maps

        if axes_correspondance == (1, 2, 3):
            ms = vg.MolSystem.from_box_data(
                resolution = res.copy(), origin = origin.copy(), deltas = vsize.copy()
            )
            obj = vg.Grid(ms, init_grid = False)
            obj.grid = data.transpose(2,1,0)
            return obj

        if axes_correspondance == (3, 2, 1):
            ms = vg.MolSystem.from_box_data(
                resolution = res[::-1], origin = origin[::-1], deltas = vsize[::-1]
            )
            obj = vg.Grid(ms, init_grid = False)
            obj.grid = data
            return obj

        raise NotImplementedError(
            f"Unsupported axes correspondence in MRC file: {axes_correspondance}. "
            "Expected (1, 2, 3) or (3, 2, 1)."
        )


# ------------------------------------------------------------------------------
