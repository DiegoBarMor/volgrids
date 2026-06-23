import numpy as np
from pathlib import Path

import volgrids as vg
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class GridIO:
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN I/O OPERATIONS
    @staticmethod
    def read_dx(path_dx: str|Path) -> "vg.Grid":
        import gridData as gd

        parser = gd.Grid(path_dx)
        box = vg.Box(parser.origin, parser.grid.shape, parser.delta)
        vgrid = vg.Grid(box, init_grid = False)
        vgrid.arr = parser.grid
        vgrid.fmt = vg.GridFormat.DX
        return vgrid


    # --------------------------------------------------------------------------
    @staticmethod
    def read_bin(path_bin: str|Path) -> "vg.Grid":
        """Load grid from a binary file in the following order (little-endian):
        - 3 unsigned ints: resolution `(nx, ny, nz)`
        - 3 floats: `deltas`
        - 3 floats: `origin`
        - `nx*ny*nz` float32 values: flattened array in C-order
        """
        import struct

        with open(path_bin, "rb") as file:
            res_bytes = file.read(12)  # 3 * 4 bytes (unsigned int)
            if len(res_bytes) < 12:
                raise ValueError("File too small to contain resolution")

            resolution = struct.unpack("<3I", res_bytes)
            deltas = struct.unpack("<3f", file.read(12))
            origin = struct.unpack("<3f", file.read(12))
            npoints = int(resolution[0]) * int(resolution[1]) * int(resolution[2])
            data = file.read(npoints * 4)
            if len(data) < npoints * 4:
                raise ValueError("File too small to contain grid data")

            box = vg.Box(origin, resolution, deltas)
            vgrid = vg.Grid(box, init_grid = False)
            vgrid.arr = np.frombuffer(data, dtype = np.float32, count = npoints).reshape(resolution)
            vgrid.fmt = vg.GridFormat.BIN
        return vgrid


    # --------------------------------------------------------------------------
    @staticmethod
    def read_mrc(path_mrc: str|Path) -> "vg.Grid":
        import gridData as gd

        with gd.mrc.mrcfile.open(path_mrc) as parser:
            ##### assume that MRC always follows the origin follows the "real space" MRC convention
            orig = parser.header["origin"]
            used_origin = np.array([orig['x'], orig['y'], orig['z']])

        vgrid = GridIO._read_mrc_ccp4(path_mrc, used_origin)
        vgrid.fmt = vg.GridFormat.MRC
        return vgrid


    # --------------------------------------------------------------------------
    @staticmethod
    def read_ccp4(path_ccp4: str|Path) -> "vg.Grid":
        import gridData as gd

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

        vgrid = GridIO._read_mrc_ccp4(path_ccp4, used_origin)
        vgrid.fmt = vg.GridFormat.CCP4
        return vgrid


    # --------------------------------------------------------------------------
    @staticmethod
    def read_cmap(path_cmap: str|Path, key: str) -> "vg.Grid":
        """Asserts that the specified key exists in the CMAP file and then reads its corresponding grid."""
        import h5py

        with h5py.File(path_cmap, 'r') as parser:
            if key not in parser["Chimera"].keys(): raise KeyError(
                f"Key '{key}' not found in '{path_cmap}'. Available keys: {list(parser['Chimera'].keys())}"
            )
            frame = parser["Chimera"][key]
            box = vg.Box(
                origin = frame.attrs["origin"],
                resolution = frame["data_zyx"].shape[::-1],
                deltas = frame.attrs["step"],
            )
            vgrid = vg.Grid(box, init_grid = False)
            vgrid.arr = frame["data_zyx"][()].transpose(2,1,0)
            vgrid.fmt = vg.GridFormat.CMAP
            vgrid.name = key
        return vgrid


    # --------------------------------------------------------------------------
    @classmethod
    def write_dx(cls, path_dx: str|Path, data: "vg.Grid"):
        path_dx = Path(path_dx)
        cls.confirm_overwrite(path_dx)
        path_dx.parent.mkdir(parents = True, exist_ok = True)

        ints = (int, np.int8, np.int16, np.int32, np.int64)
        floats = (float, np.float16, np.float32, np.float64)

        if data.arr.dtype in floats:
            grid_data = data.arr
            dtype = '"float"'
            fmt = "%.3f"
        elif data.arr.dtype in ints:
            grid_data = data.arr
            dtype = '"int"'
            fmt = "%i"
        elif data.arr.dtype == bool:
            grid_data = data.arr.astype(int)
            dtype = '"int"'
            fmt = "%i"
        else:
            raise TypeError(f"Unsupported data type for DX output: {data.arr.dtype}")

        header = '\n'.join((
            "# OpenDX density file written by volgrids",
            "# File format: http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF",
            "# Data are embedded in the header and tied to the grid positions.",
            "# Data is written in C array order: In grid[x,y,z] the axis z is fastest",
            "# varying, then y, then finally x, i.e. z is the innermost loop.",
            f"object 1 class gridpositions counts {data.xres()} {data.yres()} {data.zres()}",
            f"origin {data.xmin():6e} {data.ymin():6e} {data.zmin():6e}",
            f"delta {data.dx():6e} {0:6e} {0:6e}",
            f"delta {0:6e} {data.dy():6e} {0:6e}",
            f"delta {0:6e} {0:6e} {data.dz():6e}",
            f"object 2 class gridconnections counts  {data.xres()} {data.yres()} {data.zres()}",
            f"object 3 class array type {dtype} rank 0 items {data.npoints()}, data follows",
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
    @classmethod
    def write_bin(cls, path_bin: str|Path, data: "vg.Grid"):
        """Save grid to a binary file in the following order (little-endian):
        - 3 unsigned ints: resolution (nx, ny, nz)
        - 3 floats: deltas
        - 3 floats: origin
        - nx*ny*nz float32 values: flattened array in C-order
        """
        import struct
        path_bin = Path(path_bin)
        cls.confirm_overwrite(path_bin)
        path_bin.parent.mkdir(parents = True, exist_ok = True)

        with open(path_bin, "wb") as file:
            file.write(struct.pack("<3I", data.xres(), data.yres(), data.zres()))
            file.write(struct.pack("<3f", data.dx()  , data.dy()  , data.dz()  ))
            file.write(struct.pack("<3f", data.xmin(), data.ymin(), data.zmin()))
            file.write(data.arr.astype(np.float32).ravel(order="C").tobytes())


    # --------------------------------------------------------------------------
    @classmethod
    def write_mrc(cls, path_mrc: str|Path, data: "vg.Grid"):
        import gridData as gd

        path_mrc = Path(path_mrc)
        cls.confirm_overwrite(path_mrc)
        path_mrc.parent.mkdir(parents = True, exist_ok = True)

        with gd.mrc.mrcfile.new(path_mrc, overwrite = True) as parser:
            parser.set_data(data.arr.astype(vg.FLOAT_DTYPE).transpose(2,1,0))
            parser.voxel_size = [data.dx(), data.dy(), data.dz()]
            parser.header["origin"]['x'] = data.xmin() # MRC convention
            parser.header["origin"]['y'] = data.ymin()
            parser.header["origin"]['z'] = data.zmin()
            parser.update_header_from_data()
            parser.update_header_stats()


    # --------------------------------------------------------------------------
    @classmethod
    def write_ccp4(cls, path_ccp4: str|Path, data: "vg.Grid"):
        import gridData as gd

        path_ccp4 = Path(path_ccp4)
        cls.confirm_overwrite(path_ccp4)
        path_ccp4.parent.mkdir(parents = True, exist_ok = True)

        with gd.mrc.mrcfile.new(path_ccp4, overwrite = True) as parser:
            parser.set_data(data.arr.astype(vg.FLOAT_DTYPE).transpose(2,1,0))
            parser.voxel_size = [data.dx(), data.dy(), data.dz()]
            parser.header["origin"]['x'] = data.xmin() # MRC convention
            parser.header["origin"]['y'] = data.ymin()
            parser.header["origin"]['z'] = data.zmin()
            parser.header["nxstart"] = int(data.xmin() / data.dx()) # CCP4 convention
            parser.header["nystart"] = int(data.ymin() / data.dy())
            parser.header["nzstart"] = int(data.zmin() / data.dz())
            parser.update_header_from_data()
            parser.update_header_stats()


    # --------------------------------------------------------------------------
    @classmethod
    def write_cmap(cls, path_cmap: str|Path, data: "vg.Grid", key: str) -> None:
        import h5py
        import contextlib

        ### imitate the Chimera cmap format, as "specified" in this sample:
        ### https://github.com/RBVI/ChimeraX/blob/develop/testdata/cell15_timeseries.cmap
        def _add_generic_attrs(group, c = "GROUP"):
            group.attrs["CLASS"] = np.bytes_(c)
            group.attrs["TITLE"] = np.bytes_("")
            group.attrs["VERSION"] = np.bytes_("1.0")

        path_cmap = Path(path_cmap)
        cls.confirm_overwrite(path_cmap)
        path_cmap.parent.mkdir(parents = True, exist_ok = True)

        ### MP_CMAP_LOCK serializes concurrent writes when running trajectory MP
        lock_ctx = contextlib.nullcontext() if vg.MP_CMAP_LOCK is None else vg.MP_CMAP_LOCK
        with lock_ctx:
            if not path_cmap.exists():
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
                    frame.attrs["origin"] = data.box.min_coords.astype(vg.FLOAT_DTYPE)
                    frame.attrs["step"] = data.box.deltas.astype(vg.FLOAT_DTYPE)
                    _add_generic_attrs(frame)

                framedata = frame.create_dataset(
                    "data_zyx", data = data.arr.transpose(2,1,0), dtype = vg.FLOAT_DTYPE,
                    compression = "gzip", compression_opts = vg.GZIP_COMPRESSION
                )
                _add_generic_attrs(framedata, "CARRAY")


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ OTHER I/O UTILITIES
    @staticmethod
    def detect_format(path_grid: str|Path) -> "vg.GridFormat":
        # [TODO] improve the format detection?
        return vg.GridFormat.from_str(
            Path(path_grid).suffix.lower().lstrip('.')
        )


    # --------------------------------------------------------------------------
    @staticmethod
    def get_cmap_keys(path_cmap: Path, assert_has_keys: bool = False) -> list[str]:
        """Returns the list of keys (frame names) in a CMAP file.
        If assert_has_keys is True, raises an error if no keys are found."""
        if not path_cmap.is_file():
            if assert_has_keys:
                raise FileNotFoundError(f"CMAP file not found: {path_cmap}")
            return []

        import h5py
        with h5py.File(path_cmap, 'r') as h5:
            keys = list(h5["Chimera"].keys())

        if assert_has_keys and not keys:
            raise ValueError(f"Empty cmap file: {path_cmap}")
        return keys


    # --------------------------------------------------------------------------
    @staticmethod
    def remove(path_out: Path) -> None:
        """Removes the specified file (if it exists)."""
        path_out.unlink(missing_ok = True)


    # --------------------------------------------------------------------------
    @staticmethod
    def confirm_overwrite(path_out: Path) -> None:
        if not path_out.exists(): return
        if vg.OVERWRITE_OK: return

        response = input(f"\n\nFile {path_out} {fy.Color.red('already exists')}. Do you want to overwrite it? (y/n): ")
        if response.strip().lower() in ('y', "yes"): return

        print(fy.Color.red(f"Aborting write operation for {path_out}."))
        exit(-1)


    # --------------------------------------------------------------------------
    @staticmethod
    def restore_boolean_dtype(grid: "vg.Grid") -> "vg.Grid":
        """Restores the boolean data type of a grid that was saved as int (0 and 1) due to format limitations."""
        if not set(np.unique(grid.arr)).issubset({0, 1}): return
        grid.arr = grid.arr.astype(bool)


    # ------------------------------------------------------------------------------
    @staticmethod
    def _read_mrc_ccp4(path_mrc: Path, origin: np.ndarray) -> "vg.Grid":
        import gridData as gd

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
                box = vg.Box(origin, res, vsize)
                obj = vg.Grid(box, init_grid = False)
                obj.arr = data.transpose(2,1,0)
                return obj

            if axes_correspondance == (3, 2, 1):
                box = vg.Box(origin[::-1], res[::-1], vsize[::-1])
                obj = vg.Grid(box, init_grid = False)
                obj.arr = data
                return obj

            raise NotImplementedError(
                f"Unsupported axes correspondence in MRC file: {axes_correspondance}. "
                "Expected (1, 2, 3) or (3, 2, 1)."
            )


# //////////////////////////////////////////////////////////////////////////////
