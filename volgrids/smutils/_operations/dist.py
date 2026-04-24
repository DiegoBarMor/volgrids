import numpy as np
from pathlib import Path

import volgrids as vg

try: import freyacli as fy
except ImportError: from volgrids._vendors import freyacli as fy

PERCENTILES = (50, 75, 90, 95, 99, 99.9)

# //////////////////////////////////////////////////////////////////////////////
class DistPlot:
    @staticmethod
    def plot(path_in: Path, path_out: Path | None = None, key: str | None = None) -> None:
        """Plot non-zero voxel value distribution for a grid or CMAP trajectory frame."""
        import matplotlib
        import matplotlib.pyplot as plt

        vals, title = DistPlot._collect_values(path_in, key)

        if vals.size == 0:
            print(f"...>>> {fy.Color.red('No non-zero voxels found')} in '{path_in}'.")
            return

        print(f"Non-zero voxels: {fy.Color.yellow(str(vals.size))}")
        for p in PERCENTILES:
            # More for debugging purpose you can delete that if it bother you
            print(f"  p{p:>4}: {np.percentile(vals, p):.4g}")

        # Colorblind friendly here :) (useless with only one color but i don't care)
        matplotlib.style.use("seaborn-v0_8-colorblind")
        fig, ax = plt.subplots(figsize=(6, 3))
        # Do i add a setting to change the number of bins ? 
        ax.hist(vals, bins=80, log=True)
        ax.set_xlabel("voxel value")
        ax.set_ylabel("count (log)")
        ax.set_title(title)
        fig.tight_layout()

        if path_out is not None:
            fig.savefig(path_out, dpi=150)
            print(f"...>>> Plot saved to '{fy.Color.blue(path_out)}'.")
        else:
            plt.show()

        plt.close(fig)


    # --------------------------------------------------------------------------
    @staticmethod
    def plot_3d(
        path_in: Path, path_out: Path | None = None,
        key: str | None = None, stride: int = 2,
    ) -> None:
        """Interactive 3D volume rendering of non-zero voxels saved as HTML."""
        import plotly.graph_objects as go

        grid, title = DistPlot._collect_grid(path_in, key)
        # stride reduces resolution uniformly; the full regular grid (zeros included)
        # must be preserved so go.Volume can compute isosurfaces correctly.
        arr = grid.arr[::stride, ::stride, ::stride]
        box = grid.box

        res = arr.shape
        xs = np.linspace(box.min_coords[0], box.max_coords[0], res[0])
        ys = np.linspace(box.min_coords[1], box.max_coords[1], res[1])
        zs = np.linspace(box.min_coords[2], box.max_coords[2], res[2])
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")

        val = arr.flatten()
        nonzero_vals = val[val > 0]

        if nonzero_vals.size == 0:
            print(f"...>>> {fy.Color.red('No non-zero voxels found')} in '{path_in}'.")
            return

        print(f"Grid points: {fy.Color.yellow(str(val.size))} (stride={stride})")

        fig = go.Figure(data=go.Volume(
            x=X.flatten(), y=Y.flatten(), z=Z.flatten(),
            value=val,
            isomin=float(nonzero_vals.min()),
            isomax=float(nonzero_vals.max()),
            opacity=0.1,
            surface_count=25,
            colorscale="Plasma",
        ))
        fig.update_layout(
            title=title,
            scene_xaxis_showticklabels=False,
            scene_yaxis_showticklabels=False,
            scene_zaxis_showticklabels=False,
        )

        if path_out is not None:
            fig.write_html(str(path_out))
            print(f"...>>> 3D plot saved to '{fy.Color.blue(path_out)}'.")
        else:
            fig.show()


    # --------------------------------------------------------------------------
    @staticmethod
    def plot_echarts(
        path_in: Path, path_out: Path | None = None,
        key: str | None = None, stride: int = 1,
    ) -> None:
        """Interactive 3D scatter plot of non-zero voxels using ECharts, saved as HTML."""
        from pyecharts import options as opts
        from pyecharts.charts import Scatter3D

        grid, title = DistPlot._collect_grid(path_in, key)
        arr = grid.arr[::stride, ::stride, ::stride]
        box = grid.box

        res = arr.shape
        xs = np.linspace(box.min_coords[0], box.max_coords[0], res[0])
        ys = np.linspace(box.min_coords[1], box.max_coords[1], res[1])
        zs = np.linspace(box.min_coords[2], box.max_coords[2], res[2])
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")

        val = arr.flatten()
        x   = X.flatten()
        y   = Y.flatten()
        z   = Z.flatten()

        # Scatter3D is a point cloud — zeros can be filtered without breaking rendering
        mask = val > 0
        val, x, y, z = val[mask], x[mask], y[mask], z[mask]

        if val.size == 0:
            print(f"...>>> {fy.Color.red('No non-zero voxels found')} in '{path_in}'.")
            return

        print(f"Non-zero voxels: {fy.Color.yellow(str(val.size))} (stride={stride})")

        data = np.column_stack([x, y, z, val]).tolist()

        # plasma-like gradient (dark purple → orange → yellow)
        plasma = [
            "#0d0887", "#5302a3", "#8b0aa5", "#b83289",
            "#db5c68", "#f48849", "#febc2a", "#f0f921",
        ]

        axis_opts = opts.Axis3DOpts(
            type_="value",
            axislabel_opts=opts.LabelOpts(is_show=False),
            axisline_opts=opts.AxisLineOpts(is_show=True),
        )

        chart = (
            Scatter3D()
            .add(
                series_name="",
                data=data,
                xaxis3d_opts=axis_opts,
                yaxis3d_opts=axis_opts,
                zaxis3d_opts=axis_opts,
                grid3d_opts=opts.Grid3DOpts(width=100, height=100, depth=100),
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(title=title),
                visualmap_opts=opts.VisualMapOpts(
                    min_=float(val.min()),
                    max_=float(val.max()),
                    range_color=plasma,
                    dimension=3,
                    is_calculable=True,
                    pos_bottom="5%",
                ),
            )
        )

        if path_out is not None:
            chart.render(str(path_out))
            print(f"...>>> ECharts plot saved to '{fy.Color.blue(path_out)}'.")
        else:
            import tempfile
            import webbrowser
            with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f:
                chart.render(f.name)
                webbrowser.open(f.name)


    # --------------------------------------------------------------------------
    @staticmethod
    def _collect_grid(path_in: Path, key: str | None) -> tuple["vg.Grid", str]:
        grid = vg.GridIO.read_auto(path_in)

        if not grid.fmt.is_cmap():
            return grid, f"3D volume — {path_in.stem}"

        cmap_keys = vg.GridIO.get_cmap_keys(path_in)

        if key is not None:
            if key not in cmap_keys:
                raise ValueError(
                    f"Key '{key}' not found in '{path_in}'. "
                    f"Available keys: {cmap_keys}"
                )
            return vg.GridIO.read_cmap(path_in, key), f"3D volume — {key}"

        # average all frames into one, resampling mismatched frames to the first frame's box
        g0  = vg.GridIO.read_cmap(path_in, cmap_keys[0])
        avg = np.zeros_like(g0.arr)
        for k in cmap_keys:
            g = vg.GridIO.read_cmap(path_in, k)
            if g.box != g0.box:
                g.reshape_as_box(g0.box)
            avg += g.arr
        avg /= len(cmap_keys)

        grid_avg = vg.Grid(g0.box, init_grid=False)
        grid_avg.arr = avg
        return grid_avg, f"3D volume — {path_in.stem} (avg of {len(cmap_keys)} frames)"


    # --------------------------------------------------------------------------
    @staticmethod
    def _collect_values(path_in: Path, key: str | None) -> tuple[np.ndarray, str]:
        grid = vg.GridIO.read_auto(path_in)

        if not grid.fmt.is_cmap():
            vals = grid.arr[grid.arr > 0].ravel()
            return vals, f"Non-zero voxel distribution — {path_in.stem}"

        cmap_keys = vg.GridIO.get_cmap_keys(path_in)

        # Check if keys are alright then read
        if key is not None:
            if key not in cmap_keys:
                raise ValueError(
                    f"Key '{key}' not found in '{path_in}'. "
                    f"Available keys: {cmap_keys}"
                )
            g = vg.GridIO.read_cmap(path_in, key)
            vals = g.arr[g.arr > 0].ravel()
            # Maybe change the title, i'm not good for titles :(
            return vals, f"Non-zero voxel distribution in {key}"

        # concatenate all frames, i've not checked if it works in all weird cases
        parts = []
        for k in cmap_keys:
            g = vg.GridIO.read_cmap(path_in, k)
            parts.append(g.arr[g.arr > 0].ravel())
        vals = np.concatenate(parts) if parts else np.array([])
        # Not satisfied with this one either but idk
        title = f"Non-zero voxel distribution in {path_in.stem} ({len(cmap_keys)} frames)"
        return vals, title


# //////////////////////////////////////////////////////////////////////////////
