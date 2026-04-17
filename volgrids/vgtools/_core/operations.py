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
                resolution = f"{grid.xres()} {grid.yres()} {grid.zres()}"

            new_res = f"{grid.xres()} {grid.yres()} {grid.zres()}"
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
            if resolution is None:
                resolution = grid.box.resolution

            grid.reshape(grid.box.min_coords, grid.box.max_coords, resolution)
            vg.GridIO.write_cmap(path_out, grid, key)


    # --------------------------------------------------------------------------
    @staticmethod
    def average(path_in: Path, path_out: Path):
        keys = vg.GridIO.get_cmap_keys(path_in)
        nframes = len(keys)

        grid = vg.GridIO.read_cmap(path_in, keys[0])
        avg = np.zeros_like(grid.arr)
        for key in keys:
            print(key)
            avg += vg.GridIO.read_cmap(path_in, key).arr
        avg /= nframes

        grid_avg: vg.Grid = vg.Grid(grid.box, init_grid = False)
        grid_avg.arr = avg

        vg.GridIO.write_auto(path_out, grid_avg)


    # --------------------------------------------------------------------------
    @staticmethod
    def summary(path_in: Path):
        def numerics(g: vg.Grid, key: str):
            n_total = g.arr.size
            n_nonzero = len(g.arr[g.arr != 0])
            print(f"... grid: {key}")
            print(f"...... min: {g.arr.min():2.2e}; max: {g.arr.max():2.2e}; mean: {g.arr.mean():2.2e}")
            print(f"...... non-zero points: {n_nonzero}/{n_total} ({100*n_nonzero/n_total:.2f}%)")

        grid = vg.GridIO.read_auto(path_in)
        grid_names = vg.GridIO.get_cmap_keys(path_in) if grid.fmt.is_cmap() else [path_in.stem]

        print(f"... fmt: {grid.fmt}, ngrids: {len(grid_names)}")
        print(f"... resolution: {grid.xres()}x{grid.yres()}x{grid.zres()}; deltas: ({grid.dx():.2f},{grid.dy():.2f},{grid.dz():.2f})")
        print(f"... box: ({grid.xmin():.2f},{grid.ymin():.2f},{grid.zmin():.2f})->({grid.xmax():.2f},{grid.ymax():.2f},{grid.zmax():.2f})")

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

        deltas_0     = grid_0.box.deltas;     deltas_1     = grid_1.box.deltas
        resolution_0 = grid_0.box.resolution; resolution_1 = grid_1.box.resolution
        min_coords_0 = grid_0.box.min_coords; min_coords_1 = grid_1.box.min_coords
        max_coords_0 = grid_0.box.max_coords; max_coords_1 = grid_1.box.max_coords

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
        mask = diff.arr > threshold

        npoints_diff  = len(mask[mask])
        npoints_total = grid_0.npoints()
        cumulative_diff = np.sum(diff.arr[mask])
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
        grid.arr = vg.Math.rotate_3d(grid.arr, rotate_xy, rotate_yz, rotate_xz, in_degrees)
        vg.GridIO.write_auto(path_out, grid)


    # --------------------------------------------------------------------------
    @staticmethod
    def overlap(path_grid1: Path, path_grid2: Path, path_out: Path, operation: str = "multiply"):
        """
        Compute overlap between two molecular interaction fields.

        Args:
            path_grid1: Path to first grid file (will be interpolated)
            path_grid2: Path to second grid file (defines target coordinate system)
            path_out: Path for output grid
            operation: Type of operation ("multiply", "subtract", "add")
        """
        print(">>> Opening grids...")
        grid1 = vg.GridIO.read_auto(path_grid1)
        grid2 = vg.GridIO.read_auto(path_grid2)

        print(">>> Interpolating grid1 to match grid2 coordinate system...")
        grid1_resampled = grid1.copy()
        grid1_resampled.reshape_as(grid2)

        print(f">>> Computing overlap ({operation})...")
        result_grid = vg.Grid(grid2.ms, init_grid=False)

        if operation == "multiply":
            result_grid.arr = grid1_resampled.arr * grid2.arr
        elif operation == "subtract":
            result_grid.arr = grid1_resampled.arr - grid2.arr
        elif operation == "add":
            result_grid.arr = grid1_resampled.arr + grid2.arr
        else:
            raise ValueError(f"Unknown operation: {operation}")

        print(">>> Saving result grid...")
        path_out.parent.mkdir(parents=True, exist_ok=True)
        vg.GridIO.write_auto(path_out, result_grid)

        print(f">>> {operation.capitalize()} operation completed!")


    # --------------------------------------------------------------------------
    @staticmethod
    def overlap_cross_comparison(path_grid1: Path, path_grid2: Path, path_out: Path):
        """Cross-comparison overlap analysis (multiplication)."""
        VGOperations.overlap(path_grid1, path_grid2, path_out, "multiply")


    # --------------------------------------------------------------------------
    @staticmethod
    def overlap_difference(path_grid1: Path, path_grid2: Path, path_out: Path):
        """Difference overlap analysis (subtraction)."""
        VGOperations.overlap(path_grid1, path_grid2, path_out, "subtract")


# //////////////////////////////////////////////////////////////////////////////
