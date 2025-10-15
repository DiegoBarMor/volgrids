import numpy as np
from pathlib import Path

import volgrids as vg
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
                resolution = f"{grid.xres} {grid.yres} {grid.zres}"

            new_res = f"{grid.xres} {grid.yres} {grid.zres}"
            if (new_res != resolution) and not warned:
                print(
                    f">>> Warning: Grid {path_in} has different resolution {new_res} than the first grid {resolution}. " +\
                    "Chimera won't recognize it as a volume series and open every grid in a separate representation." +\
                    "Use `run/vgtools.py fix_cmap` if you want to fix this."
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


    # --------------------------------------------------------------------------
    @staticmethod
    def compare(path_in_0: Path, path_in_1: Path, threshold: float) -> "vgt.ComparisonResult":
        def _are_different_vector(vec0, vec1):
            diff = np.abs(vec0 - vec1)
            return len(diff[diff > threshold]) != 0

        grid_0 = vg.GridIO.read_auto(path_in_0)
        grid_1 = vg.GridIO.read_auto(path_in_1)

        deltas_0     = grid_0.get_deltas();     deltas_1     = grid_1.get_deltas()
        resolution_0 = grid_0.get_resolution(); resolution_1 = grid_1.get_resolution()
        min_coords_0 = grid_0.get_min_coords(); min_coords_1 = grid_1.get_min_coords()
        max_coords_0 = grid_0.get_max_coords(); max_coords_1 = grid_1.get_max_coords()

        if _are_different_vector(resolution_0, resolution_1):
            return vgt.ComparisonResult(0, 0, 0.0, 0.0,
                [f"Warning: Grids {path_in_0} and {path_in_1} have different shapes: {resolution_0} vs {resolution_1}. Aborting."]
            )

        warnings = []
        if _are_different_vector(min_coords_0, min_coords_1):
            warnings.append(
                f"Warning: Grids {path_in_0} and {path_in_1} have different origin: {min_coords_0} vs {min_coords_1}. Comparison may not be accurate."
            )
        if _are_different_vector(max_coords_0, max_coords_1):
            warnings.append(
                f"Warning: Grids {path_in_0} and {path_in_1} have different max coordinate: {max_coords_0} vs {max_coords_1}. Comparison may not be accurate."
            )
        if _are_different_vector(deltas_0, deltas_1):
            warnings.append(
                f"Warning: Grids {path_in_0} and {path_in_1} have different deltas: {deltas_0} vs {deltas_1}. Comparison may not be accurate."
            )

        diff = abs(grid_1 - grid_0)
        mask = diff.grid > threshold

        npoints_diff  = len(mask[mask])
        npoints_total = grid_0.xres * grid_0.yres * grid_0.zres
        cumulative_diff = np.sum(diff.grid[mask])
        avg_diff = (cumulative_diff / npoints_diff) if (npoints_diff > 0) else 0

        return vgt.ComparisonResult(npoints_diff, npoints_total, cumulative_diff, avg_diff, warnings)



# //////////////////////////////////////////////////////////////////////////////
