import sys, os, h5py
import numpy as np
import gridData as gd
from pathlib import Path
import volgrids.vgrids as vg

################################################################################
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN I/O OPERATIONS
def read_mrc(path_mrc) -> "vg.Grid":
    with gd.mrc.mrcfile.open(path_mrc) as grid_mrc:
        vsize = np.array(grid_mrc.voxel_size)
        origin = np.array(grid_mrc.header["origin"])

        grid_vg = vg.Grid(**_grid_init_metadata(
            resolution = np.array(grid_mrc.data.shape),
            origin = np.array([origin["x"], origin["y"], origin["z"]], dtype = vg.FLOAT_DTYPE),
            delta = np.array([vsize["x"], vsize["y"], vsize["z"]], dtype = vg.FLOAT_DTYPE)
        ))
        grid_vg.grid = grid_mrc.data.transpose(2,1,0)
    return grid_vg


# ------------------------------------------------------------------------------
def read_dx(path_dx) -> "vg.Grid":
    grid_dx = gd.Grid(path_dx)
    grid_vg = vg.Grid(**_grid_init_metadata(
        resolution = grid_dx.grid.shape,
        origin = grid_dx.origin,
        delta = grid_dx.delta
    ))
    grid_vg.grid = grid_dx.grid
    return grid_vg


# ------------------------------------------------------------------------------
def read_cmap(path_cmap, key) -> "vg.Grid":
    with h5py.File(path_cmap, 'r') as h5:
        frame = h5["Chimera"][key]
        rz, ry, rx = frame["data_zyx"].shape
        ox, oy, oz = frame.attrs["origin"]
        dz, dy, dx = frame.attrs["step"]
        grid_vg = vg.Grid(**_grid_init_metadata(
            resolution = np.array([rx, ry, rz]),
            origin = np.array([ox, oy, oz]),
            delta = np.array([dx, dy, dz])
        ))
        grid_vg.grid = frame["data_zyx"][()].transpose(2,1,0)
    return grid_vg


# ------------------------------------------------------------------------------
def write_mrc(path_mrc, data: "vg.Grid"):
    with gd.mrc.mrcfile.new(path_mrc, overwrite = True) as grid_mrc:
        grid_mrc.set_data(data.grid.transpose(2,1,0))
        grid_mrc.voxel_size = [data.dx, data.dy, data.dz]
        grid_mrc.header["origin"]['x'] = data.xmin
        grid_mrc.header["origin"]['y'] = data.ymin
        grid_mrc.header["origin"]['z'] = data.zmin
        grid_mrc.update_header_from_data()
        grid_mrc.update_header_stats()


# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
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


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ OTHER I/O UTILITIES
def read_auto(path_grid: Path) -> tuple[str, "vg.Grid"]:
    """Detect the format of the grid file based on its extension and then read it."""
    ext = path_grid.suffix.lower()

    if ext == ".dx":
        return "DX", read_dx(path_grid)

    if ext == ".mrc":
        return "MRC", read_mrc(path_grid)

    elif ext == ".cmap":
        keys = get_cmap_keys(path_grid)
        if not keys: raise ValueError(f"Empty cmap file: {path_grid}")
        return "CMAP", read_cmap(path_grid, keys[0])

    raise ValueError(f"Unrecognized file format: {ext}")


# ------------------------------------------------------------------------------
def get_cmap_keys(path_cmap) -> list[str]:
    with h5py.File(path_cmap, 'r') as h5:
        return list(h5["Chimera"].keys())


# ------------------------------------------------------------------------------
def resolve_path(path: Path):
    """Resolve the path to the project root directory."""
    project_root = Path(__file__).resolve().parent.parent.parent
    return project_root / path


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


################################################################################
