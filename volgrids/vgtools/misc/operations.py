from pathlib import Path

import volgrids.vgrids as vg
import volgrids.vgtools as vgt

# //////////////////////////////////////////////////////////////////////////////
class VGOperations:
    @staticmethod
    def convert(path_in: Path, path_out: Path, fmt_out: vg.GridFormat):
        grid = vg.GridIO.read_auto(vgt.PATH_CONVERT_IN)

        func: callable = {
            vg.GridFormat.DX: vg.GridIO.write_dx,
            vg.GridFormat.MRC: vg.GridIO.write_mrc,
            vg.GridFormat.CCP4: vg.GridIO.write_ccp4,
            vg.GridFormat.CMAP: vg.GridIO.write_cmap,
        }.get(fmt_out, None)
        if func is None:
            raise ValueError(f"Unknown format for conversion: {fmt_out}")

        extra_args = (path_in.stem,) if fmt_out == vg.GridFormat.CMAP else ()
        func(path_out, grid, *extra_args)


    # --------------------------------------------------------------------------
    @staticmethod
    def pack(paths_in: list[Path], path_out: Path):
        resolution = None
        warned = False
        for path_in in paths_in:
            grid = vg.GridIO.read_auto(path_in)
            if resolution is None:
                resolution = (grid.xres, grid.yres, grid.zres)

            new_res = (grid.xres, grid.yres, grid.zres)
            if (new_res != resolution) and not warned:
                print(
                    f">>> Warning: Grid {path_in} has different resolution {new_res} than the first grid {resolution}. " +\
                    "Chimera won't recognize it as a volume series and open every grid in a separate representation." +\
                    "Use `smiffer.py fix_cmap` if you want to fix this."
                )
                warned = True

            key = str(path_in.parent / path_in.stem).replace(' ', '_').replace('/', '_').replace('\\', '_')
            # key = path_in.stem
            vg.GridIO.write_cmap(path_out, grid, key)


    # --------------------------------------------------------------------------
    @staticmethod
    def unpack(path_in: Path, folder_out: Path):
        keys = vg.GridIO.get_cmap_keys(path_in)
        for key in keys:
            path_out = folder_out / f"{key}.cmap"
            grid = vg.GridIO.read_cmap(path_in, key)
            vg.GridIO.write_cmap(path_out, grid, key)


    # --------------------------------------------------------------------------
    @staticmethod
    def fix_cmap(path_in: Path, path_out: Path):
        resolution = None
        keys = vg.GridIO.get_cmap_keys(path_in)
        for key in keys:
            grid = vg.GridIO.read_cmap(path_in, key)

            minCoords = (grid.xmin, grid.ymin, grid.zmin)
            maxCoords = (grid.xmax, grid.ymax, grid.zmax)
            if resolution is None:
                resolution = (grid.xres, grid.yres, grid.zres)

            grid.reshape(minCoords, maxCoords, resolution)
            vg.GridIO.write_cmap(path_out, grid, key)


# //////////////////////////////////////////////////////////////////////////////
