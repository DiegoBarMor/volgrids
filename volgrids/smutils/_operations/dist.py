import numpy as np
from pathlib import Path

import volgrids as vg

try: import freyacli as fy
except ImportError: from volgrids._vendors import freyacli as fy


# //////////////////////////////////////////////////////////////////////////////
class Histogram:
    PERCENTILES = (50, 75, 90, 95, 99, 99.9)

    # --------------------------------------------------------------------------
    @classmethod
    def plot(cls, path_in: Path, path_out: Path | None = None, key: str | None = None) -> None:
        """Plot non-zero voxel value distribution for a grid or CMAP trajectory frame."""
        import matplotlib
        import matplotlib.pyplot as plt
        matplotlib.style.use("seaborn-v0_8-colorblind")

        vals, title = cls._collect_values(path_in, key)

        if vals.size == 0:
            print(f"...>>> {fy.Color.red('All voxels are 0')} in '{fy.Color.blue(path_in)}'.")
            return

        print(f"Non-zero voxels: {fy.Color.yellow(vals.size)}")
        ### [TODO] print as percentage too, similar to `vgtools summary`

        ### [TODO] toggling detailed print could be regulated by a config, will leave like this for now
        for p in cls.PERCENTILES:
            str_p = fy.Color.green(f"{p:>4}")
            print(f"  p{str_p}: {np.percentile(vals, p):.4g}")

        fig, ax = plt.subplots(figsize=(6, 3))
        ax.hist(vals, bins=80, log=True) # [TODO] add config for bins
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
        def get_nonzero(path_in, key = "") -> np.ndarray:
            grid = vg.GridIO.read_cmap(path_in, key) \
                if key else vg.GridIO.read_auto(path_in)
            return grid.arr[grid.arr > 0].ravel()

        fmt = vg.GridIO.detect_format(path_in)
        str_preffix = "Non-zero voxel distribution"

        if not fmt.is_cmap():
            return get_nonzero(path_in), f"{str_preffix}: {path_in.stem}"

        if key is None:
            cmap_keys = vg.GridIO.get_cmap_keys(path_in)
            parts = [get_nonzero(path_in, k) for k in cmap_keys]
            vals = np.concatenate(parts) if parts else np.array([]) # [TODO] test with large trajectores
            title = f"{str_preffix}: {path_in.stem} ({len(cmap_keys)} frames)"
            return vals, title

        return get_nonzero(path_in, key), f"{str_preffix}, key: {key}"


# //////////////////////////////////////////////////////////////////////////////
