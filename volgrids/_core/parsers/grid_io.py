import os, h5py
import numpy as np
import gridData as gd
import volgrids as vg
from pathlib import Path
from enum import Enum, auto

# //////////////////////////////////////////////////////////////////////////////
class GridFormat(Enum):
    DX = auto()
    MRC = auto()
    CCP4 = auto()
    CMAP = auto()
    CMAP_PACKED = auto()

    @classmethod
    def from_str(cls, value: str) -> "GridFormat":
        """Parse a string to a GridFormat enum."""
        k = value.upper()
        if k not in cls.__members__:
            raise ValueError(f"Invalid output format: {k}. Valid options are: {list(cls.__members__.keys())}")
        return cls[k]


# //////////////////////////////////////////////////////////////////////////////
class GridIO:
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN I/O OPERATIONS
    @staticmethod
    def read_dx(path_dx) -> "vg.Grid":
        parser_dx = gd.Grid(path_dx)
        obj = vg.Grid(**_grid_init_metadata(
            resolution = parser_dx.grid.shape,
            origin = parser_dx.origin,
            delta = parser_dx.delta
        ))
        obj.grid = parser_dx.grid
        return obj


    # --------------------------------------------------------------------------
    @staticmethod
    def read_mrc(path_mrc) -> "vg.Grid":
        with gd.mrc.mrcfile.open(path_mrc) as parser_mrc:
            # machine_stamp = parser_mrc.header.machst
            ### [68 68 0 0] or [68 65 0 0] for little-endian <--- tested
            ### [17 17 0 0] for big-endian <--- what happens in these cases?

            axes_correspondance =\
                parser_mrc.header.mapc, parser_mrc.header.mapr, parser_mrc.header.maps

            vsize = np.array([
                parser_mrc.voxel_size['x'],
                parser_mrc.voxel_size['y'],
                parser_mrc.voxel_size['z'],
            ], dtype = vg.FLOAT_DTYPE)

            origin = np.array([
                parser_mrc.header["origin"]['x'] + vsize[0] * parser_mrc.header["nxstart"],
                parser_mrc.header["origin"]['y'] + vsize[1] * parser_mrc.header["nystart"],
                parser_mrc.header["origin"]['z'] + vsize[2] * parser_mrc.header["nzstart"],
            ], dtype = vg.FLOAT_DTYPE)

            if axes_correspondance == (1, 2, 3):
                obj = vg.Grid(**_grid_init_metadata(
                    resolution = np.array([parser_mrc.header.mx, parser_mrc.header.my, parser_mrc.header.mz]),
                    origin = origin.copy(),
                    delta = vsize.copy()
                ))
                obj.grid = parser_mrc.data.transpose(2,1,0)

            elif axes_correspondance == (3, 2, 1):
                obj = vg.Grid(**_grid_init_metadata(
                    resolution = np.array([parser_mrc.header.mz, parser_mrc.header.my, parser_mrc.header.mx]),
                    origin = origin[::-1],
                    delta = vsize[::-1]
                ))
                obj.grid = parser_mrc.data.copy()

            else:
                raise NotImplementedError(
                    f"Unsupported axes correspondence in MRC file: {axes_correspondance}. "
                    "Expected (1, 2, 3) or (3, 2, 1)."
                )

        return obj


    # --------------------------------------------------------------------------
    @staticmethod
    def read_ccp4(path_ccp4) -> "vg.Grid":
        return GridIO.read_mrc(path_ccp4)


    # --------------------------------------------------------------------------
    @staticmethod
    def read_cmap(path_cmap, key) -> "vg.Grid":
        with h5py.File(path_cmap, 'r') as parser_cmap:
            frame = parser_cmap["Chimera"][key]
            rz, ry, rx = frame["data_zyx"].shape
            ox, oy, oz = frame.attrs["origin"]
            dz, dy, dx = frame.attrs["step"]
            obj = vg.Grid(**_grid_init_metadata(
                resolution = np.array([rx, ry, rz]),
                origin = np.array([ox, oy, oz]),
                delta = np.array([dx, dy, dz])
            ))
            obj.grid = frame["data_zyx"][()].transpose(2,1,0)
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
        with gd.mrc.mrcfile.new(path_mrc, overwrite = True) as grid_mrc:
            grid_mrc.set_data(data.grid.transpose(2,1,0))
            grid_mrc.voxel_size = [data.dx, data.dy, data.dz]
            grid_mrc.header["origin"]['x'] = data.xmin
            grid_mrc.header["origin"]['y'] = data.ymin
            grid_mrc.header["origin"]['z'] = data.zmin
            grid_mrc.update_header_from_data()
            grid_mrc.update_header_stats()


    # --------------------------------------------------------------------------
    @staticmethod
    def write_ccp4(path_ccp4, data: "vg.Grid"):
        GridIO.write_mrc(path_ccp4, data)


    # --------------------------------------------------------------------------
    @staticmethod
    def write_cmap(path_cmap, data: "vg.Grid", key):
        ### imitate the Chimera cmap format, as "specified" in this sample:
        ### https://github.com/RBVI/ChimeraX/blob/develop/testdata/cell15_timeseries.cmap
        def add_generic_attrs(group, c = "GROUP"):
            group.attrs["CLASS"] = np.bytes_(c)
            group.attrs["TITLE"] = np.bytes_("")
            group.attrs["VERSION"] = np.bytes_("1.0")

        if not os.path.exists(path_cmap):
            with h5py.File(path_cmap, 'w') as h5:
                h5.attrs["PYTABLES_FORMAT_VERSION"] = np.bytes_("2.0")
                add_generic_attrs(h5)

                chim = h5.create_group("Chimera")
                add_generic_attrs(chim)

        with h5py.File(path_cmap, 'a') as h5:
            chim = h5["Chimera"]
            if key in chim.keys():
                frame = chim[key]
                if "data_zyx" in frame.keys():
                    del frame["data_zyx"]
            else:
                frame = h5.create_group(f"/Chimera/{key}")
                frame.attrs["chimera_map_version"] = np.int64(1)
                frame.attrs["chimera_version"] = np.bytes_(b'1.12_b40875')
                frame.attrs["name"] = np.bytes_(key)
                frame.attrs["origin"] = np.array([data.xmin, data.ymin, data.zmin], dtype = vg.FLOAT_DTYPE)
                frame.attrs["step"] = np.array([data.dz, data.dy, data.dx], dtype = vg.FLOAT_DTYPE)
                add_generic_attrs(frame)

            framedata = frame.create_dataset(
                "data_zyx", data = data.grid.transpose(2,1,0), dtype = vg.FLOAT_DTYPE,
                compression = "gzip", compression_opts = vg.GZIP_COMPRESSION
            )
            add_generic_attrs(framedata, "CARRAY")


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ OTHER I/O UTILITIES
    @staticmethod
    def read_auto(path_grid: Path) -> tuple[GridFormat, "vg.Grid"]:
        """Detect the format of the grid file based on its extension and then read it."""
        ext = path_grid.suffix.lower()

        # [TODO] improve the format detection?
        if ext == ".dx":
            return GridFormat.DX, GridIO.read_dx(path_grid)

        if ext == ".mrc":
            return GridFormat.MRC, GridIO.read_mrc(path_grid)

        if ext == ".ccp4":
            return GridFormat.CCP4, GridIO.read_ccp4(path_grid)

        if ext == ".cmap":
            keys = GridIO.get_cmap_keys(path_grid)
            if not keys: raise ValueError(f"Empty cmap file: {path_grid}")
            return GridFormat.CMAP, GridIO.read_cmap(path_grid, keys[0])

        raise ValueError(f"Unrecognized file format: {ext}")


    # --------------------------------------------------------------------------
    @staticmethod
    def get_cmap_keys(path_cmap) -> list[str]:
        with h5py.File(path_cmap, 'r') as h5:
            return list(h5["Chimera"].keys())


# //////////////////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------------------
def _grid_init_metadata(resolution: np.array, origin: np.array, delta: np.array)-> dict:
    return dict(
        data = {
            "resolution": resolution,
            "minCoords": origin,
            "maxCoords": origin + delta * resolution,
            "deltas": delta
        },
        init_grid = False
    )


# ------------------------------------------------------------------------------
