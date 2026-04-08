import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.vgtools as vgt

# //////////////////////////////////////////////////////////////////////////////
class VGOperations:
    @staticmethod
    def convert(path_in: Path, path_out: Path, fmt_out: vg.GridFormat):
        grid = vg.GridIO.read_auto(path_in)

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
                    "Chimera won't recognize it as a volume series and open every grid in a separate representation. " +\
                    "Use `volgrids vgtools fix_cmap` if you want to fix this."
                )
                warned = True

            key = str(path_in.parent / path_in.stem).replace(' ', '_').replace('/', '_').replace('\\', '_')
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
    def average(path_in: Path, path_out: Path):
        keys = vg.GridIO.get_cmap_keys(path_in)
        nframes = len(keys)

        grid = vg.GridIO.read_cmap(path_in, keys[0])
        avg = np.zeros_like(grid.grid)
        for key in keys:
            print(key)
            avg += vg.GridIO.read_cmap(path_in, key).grid
        avg /= nframes

        grid_avg: vg.Grid = vg.Grid(grid.ms, init_grid = False)
        grid_avg.grid = avg

        vg.GridIO.write_auto(path_out, grid_avg)


    # --------------------------------------------------------------------------
    @staticmethod
    def summary(path_in: Path):
        def numerics(g: vg.Grid, key: str):
            n_total = g.grid.size
            n_nonzero = len(g.grid[g.grid != 0])
            print(f"... grid: {key}")
            print(f"...... min: {g.grid.min():2.2e}; max: {g.grid.max():2.2e}; mean: {g.grid.mean():2.2e}")
            print(f"...... non-zero points: {n_nonzero}/{n_total} ({100*n_nonzero/n_total:.2f}%)")

        grid = vg.GridIO.read_auto(path_in)
        grid_names = vg.GridIO.get_cmap_keys(path_in) if grid.fmt.is_cmap() else [path_in.stem]

        print(f"... fmt: {grid.fmt}, ngrids: {len(grid_names)}")
        print(f"... resolution: {grid.xres}x{grid.yres}x{grid.zres}; deltas: ({grid.dx:.2f},{grid.dy:.2f},{grid.dz:.2f})")
        print(f"... box: ({grid.xmin:.2f},{grid.ymin:.2f},{grid.zmin:.2f})->({grid.xmax:.2f},{grid.ymax:.2f},{grid.zmax:.2f})")

        if not grid.fmt.is_cmap():
            numerics(grid, path_in.stem); print()
            return

        for key in grid_names:
            numerics(vg.GridIO.read_cmap(path_in, key), key)
        print()


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


    # --------------------------------------------------------------------------
    @staticmethod
    def rotate(
        path_in: Path, path_out: Path,
        rotate_xy: float, rotate_yz: float, rotate_xz: float,
        in_degrees: bool = True
    ) -> None:
        grid = vg.GridIO.read_auto(path_in)
        vg.GridIO.restore_boolean_dtype(grid)
        grid.grid = vg.Math.rotate_3d(grid.grid, rotate_xy, rotate_yz, rotate_xz, in_degrees)
        vg.GridIO.write_auto(path_out, grid)


# //////////////////////////////////////////////////////////////////////////////
