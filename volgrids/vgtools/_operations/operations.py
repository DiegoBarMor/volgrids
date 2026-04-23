import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.vgtools as vgt

try: import freyacli as fy # to display colored text
except ImportError: from volgrids._vendors import freyacli as fy

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
                    f">>> {fy.Color.red('Warning')}: Grid {path_in} has different resolution {new_res} than the first grid {resolution}. " +\
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

            str_min_ = fy.Color.red(f"{g.arr.min():.2e}")
            str_max_ = fy.Color.blue(f"{g.arr.max():.2e}")
            str_perc = fy.Color.yellow(f"{100*n_nonzero/n_total:.2f}%")

            print(f"... grid: {fy.Color.cyan(key)}")
            print(f"...... min: {str_min_}; max: {str_max_}; mean: {g.arr.mean():2.2e}")
            print(f"...... non-zero points: {n_nonzero}/{n_total} ({str_perc})")

        grid = vg.GridIO.read_auto(path_in)
        grid_names = vg.GridIO.get_cmap_keys(path_in) if grid.fmt.is_cmap() else [path_in.stem]

        str_res = fy.Color.green(f"{grid.xres()}x{grid.yres()}x{grid.zres()}")
        str_min = fy.Color.red(f"{grid.xmin():.2f},{grid.ymin():.2f},{grid.zmin():.2f}")
        str_max = fy.Color.blue(f"{grid.xmax():.2f},{grid.ymax():.2f},{grid.zmax():.2f}")

        print(f"... fmt: {fy.Color.magenta(grid.fmt.name)}, ngrids: {fy.Color.yellow(len(grid_names))}")
        print(f"... resolution: {str_res}; deltas: ({grid.dx():.2f},{grid.dy():.2f},{grid.dz():.2f})")
        print(f"... box: ({str_min})->({str_max})")

        if not grid.fmt.is_cmap():
            numerics(grid, path_in.stem); print()
            return

        for key in grid_names:
            numerics(vg.GridIO.read_cmap(path_in, key), key)
        print()


    # --------------------------------------------------------------------------
    @staticmethod
    def compare(path_in_0: Path, path_in_1: Path, threshold: float) -> "vgt.ComparisonResult":
        grid_0 = vg.GridIO.read_auto(path_in_0)
        grid_1 = vg.GridIO.read_auto(path_in_1)

        deltas_0     = grid_0.box.deltas;     deltas_1     = grid_1.box.deltas
        resolution_0 = grid_0.box.resolution; resolution_1 = grid_1.box.resolution
        min_coords_0 = grid_0.box.min_coords; min_coords_1 = grid_1.box.min_coords
        max_coords_0 = grid_0.box.max_coords; max_coords_1 = grid_1.box.max_coords

        str_warning = fy.Color.red("Warning")

        if not np.allclose(resolution_0, resolution_1):
            return vgt.ComparisonResult(0, 0, 0.0, 0.0,
                [f"{str_warning}: Grids {path_in_0} and {path_in_1} have different shapes: {resolution_0} vs {resolution_1}. {fy.Color.red('Aborting', bright = False)}."]
            )

        warnings = []
        if not np.allclose(min_coords_0, min_coords_1):
            warnings.append(
                f"{str_warning}: Grids {path_in_0} and {path_in_1} have different origin: {min_coords_0} vs {min_coords_1}. Comparison may not be accurate."
            )
        if not np.allclose(max_coords_0, max_coords_1):
            warnings.append(
                f"{str_warning}: Grids {path_in_0} and {path_in_1} have different max coordinate: {max_coords_0} vs {max_coords_1}. Comparison may not be accurate."
            )
        if not np.allclose(deltas_0, deltas_1):
            warnings.append(
                f"{str_warning}: Grids {path_in_0} and {path_in_1} have different deltas: {deltas_0} vs {deltas_1}. Comparison may not be accurate."
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
    def op(
        operation: callable, path_out: Path, path_in_0: Path, path_in_1: Path = None,
    ) -> None:
        """
        Perform the numeric `operation` between two grids.
        Having `path_in_1` be None implies a unary operation (e.g. abs) over `path_in_0`.
        Supports multi-frame CMAP trajectories: frames are processed one-by-one, with
        broadcasting if one side has a single frame and the other has N.
        """
        path_in_0 = Path(path_in_0)
        path_out  = Path(path_out)

        if path_in_1 is None: # unary operation
            if path_in_0.suffix.lower() == ".cmap":
                keys = vg.GridIO.get_cmap_keys(path_in_0)
                if not keys: raise ValueError(f"Empty cmap file: {path_in_0}")
                path_out.parent.mkdir(parents=True, exist_ok=True)
                for key in keys:
                    grid = vg.GridIO.read_cmap(path_in_0, key)
                    vg.GridIO.write_cmap(path_out, operation(grid), key=key)
            else:
                grid = vg.GridIO.read_auto(path_in_0)
                vg.GridIO.write_auto(path_out, operation(grid))
            return

        path_in_1 = Path(path_in_1)
        ext0, ext1 = path_in_0.suffix.lower(), path_in_1.suffix.lower()

        if ext0 == ".cmap" and ext1 == ".cmap":
            keys0 = vg.GridIO.get_cmap_keys(path_in_0)
            keys1 = vg.GridIO.get_cmap_keys(path_in_1)
            if not keys0: raise ValueError(f"Empty cmap file: {path_in_0}")
            if not keys1: raise ValueError(f"Empty cmap file: {path_in_1}")

            n0, n1 = len(keys0), len(keys1)
            if n0 != n1 and n0 != 1 and n1 != 1:
                raise ValueError(
                    f"Incompatible trajectory lengths: {path_in_0} has {n0} frames, "
                    f"{path_in_1} has {n1} frames. Must be equal or one must be 1."
                )

            n_frames = max(n0, n1)
            path_out.parent.mkdir(parents=True, exist_ok=True)
            for i in range(n_frames):
                k0 = keys0[i if n0 > 1 else 0]
                k1 = keys1[i if n1 > 1 else 0]
                grid_0 = vg.GridIO.read_cmap(path_in_0, k0)
                grid_1 = vg.GridIO.read_cmap(path_in_1, k1)
                if grid_0.box != grid_1.box:
                    str_grid_0 = f"'{fy.Color.blue(path_in_0)}' {fy.Color.yellow(grid_0.box.resolution)}"
                    str_grid_1 = f"'{fy.Color.red(path_in_1)}' {fy.Color.yellow(grid_1.box.resolution)}"
                    print(f"...>>> Interpolating {str_grid_1} to match {str_grid_0} coordinate system...")
                    grid_1.reshape_as_box(grid_0.box)
                vg.GridIO.write_cmap(path_out, operation(grid_0, grid_1), key=k0)
            return

        grid_0 = vg.GridIO.read_auto(path_in_0)
        grid_1 = vg.GridIO.read_auto(path_in_1)

        if grid_0.box != grid_1.box:
            str_grid_0 = f"'{fy.Color.blue(path_in_0)}' {fy.Color.yellow(grid_0.box.resolution)}"
            str_grid_1 = f"'{fy.Color.red(path_in_1)}' {fy.Color.yellow(grid_1.box.resolution)}"
            print(f"...>>> Interpolating {str_grid_1} to match {str_grid_0} coordinate system...")
            grid_1.reshape_as_box(grid_0.box)

        vg.GridIO.write_auto(path_out, operation(grid_0, grid_1))


# //////////////////////////////////////////////////////////////////////////////
