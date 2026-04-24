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
