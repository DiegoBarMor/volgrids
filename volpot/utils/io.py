import os, json
import numpy as np
import volpot as vp
import gridData as gd

################################################################################
def read_json(path_json) -> dict | list:
    with open(path_json, 'r') as file:
        return json.load(file)


# ------------------------------------------------------------------------------
def read_mrc(path_mrc) -> "vp.Grid":
    with gd.mrc.mrcfile.open(path_mrc) as grid_mrc:
        vsize = np.array(grid_mrc.voxel_size)
        origin = np.array(grid_mrc.header["origin"])

        delta = np.array([vsize["x"], vsize["y"], vsize["z"]], dtype = np.float32)
        resolution = np.array(grid_mrc.data.shape)
        minCoords = np.array([origin["x"], origin["y"], origin["z"]], dtype = np.float32)

        grid_vp = vp.Grid(
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


# ------------------------------------------------------------------------------
def read_dx(path_dx) -> "vp.Grid":
    grid_dx = gd.Grid(path_dx)
    grid_vp = vp.Grid(
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


# ------------------------------------------------------------------------------
def write_json(path_json, data: dict | list):
    with open(path_json, 'w') as file:
        file.write(json.dumps(data))


# --------------------------------------------------------------------------
def write_mrc(path_mrc, data: "vp.Grid"):
    with gd.mrc.mrcfile.new(path_mrc, overwrite = True) as grid_mrc:
        grid_mrc.set_data(data.grid.T.astype(np.float32))
        grid_mrc.voxel_size = [data.dx, data.dy, data.dz]
        grid_mrc.header["origin"]['x'] = data.xmin
        grid_mrc.header["origin"]['y'] = data.ymin
        grid_mrc.header["origin"]['z'] = data.zmin
        grid_mrc.update_header_from_data()
        grid_mrc.update_header_stats()


# ------------------------------------------------------------------------------
def write_dx(path_dx, data: "vp.Grid"):
    if data.grid.dtype == np.float32:
        grid_data = data.grid
        dtype = '"float"'
        fmt = "%.3f"
    elif data.grid.dtype == int:
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
        "# OpenDX density file written by volpot2",
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
def save_metadata(metadata):
    meta = metadata.copy()
    for k in ("pdb", "out", "apbs", "meta"):
        meta[k] = str(meta[k])

    path_json = metadata["meta"]
    os.makedirs(path_json.parent, exist_ok = True)
    write_json(path_json, meta)


################################################################################
