import numpy as np
import pandas as pd

import volgrids as vg
import volgrids.veins as ve

# ------------------------------------------------------------------------------
def _assert_df_has_cols(df: pd.DataFrame, *cols_metadata):
    if not set(df.columns).issuperset(cols_metadata):
        raise ValueError(
            f"CSV file '{ve.PATH_CSV_IN}' must contain the columns: " +\
            ", ".join(map(lambda x: f"'{x}'", cols_metadata)) + " " +\
            f"Found columns: {df.columns}"
        )


# //////////////////////////////////////////////////////////////////////////////
class AppVeins(vg.App):
    CONFIG_MODULES = (vg,)
    _CLASS_PARAM_HANDLER = ve.ParamHandlerVeins
    COLS_POSITION = ["xpos", "ypos", "zpos"]

    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.df = pd.read_csv(ve.PATH_CSV_IN).dropna(how = "any")
        self.cols_values = set(self.df.columns) - set(self.COLS_POSITION)

        _assert_df_has_cols(self.df, *self.COLS_POSITION)
        self.coords = self.df[self.COLS_POSITION].values.astype(float)

        self.ms = vg.MolSystem.from_box_data(
            origin    = np.min(self.coords, axis = 0) - vg.EXTRA_BOX_SIZE,
            maxCoords = np.max(self.coords, axis = 0) + vg.EXTRA_BOX_SIZE,
            molname = ve.PATH_CSV_IN.name
        )

        self.timer = vg.Timer(
            f">>> Now processing '{ve.PATH_CSV_IN.name}' ({ve.MODE})"
        )


    # --------------------------------------------------------------------------
    def run(self):
        title_suffix = ve.MODE[:5] +\
            ('' if ve.TRAJ_FRAME is None else f"_frame{ve.TRAJ_FRAME:04}")
        process_grids: callable = self._get_process_grids_func()

        self.timer.start()
        process_grids(title_suffix)
        self.timer.end()


    # --------------------------------------------------------------------------
    def _import_config_dependencies(self):
        import numpy as np
        return {"np": np, "vg": vg}


    # --------------------------------------------------------------------------
    def _get_process_grids_func(self) -> callable:
        match ve.MODE:
            case "energy": return self._process_grids_energy
            case "force":  return self._process_grids_force
            case _: raise ValueError(f"Invalid mode '{ve.MODE}'.")


    # --------------------------------------------------------------------------
    def _process_grids_energy(self, title_suffix: str):
        raise NotImplementedError("[WIP]")
        for kind in self.df["kind"].unique():
            grid = ve.GridVolumetricEnergy(self.ms, self.df, kind) # [WIP]
            grid.populate_grid() # [WIP]
            grid.save_data(ve.FOLDER_OUT, f"{grid.kind}{title_suffix}")


    # --------------------------------------------------------------------------
    def _process_grids_force(self, title_suffix: str):
        force_kinds = self._infer_force_column_names()
        for kind in force_kinds:
            forces = self.df[[f"x{kind}", f"y{kind}", f"z{kind}"]].values.astype(float)
            grid = ve.GridVolumetricForce(self.ms, self.coords, forces)
            grid.populate_grid()
            grid.save_data(ve.FOLDER_OUT, f"{kind}{title_suffix}")


    # --------------------------------------------------------------------------
    def _infer_force_column_names(self) -> list[str]:
        cols_start_x = set(c[1:] for c in self.cols_values if c.startswith('x'))
        cols_start_y = set(c[1:] for c in self.cols_values if c.startswith('y'))
        cols_start_z = set(c[1:] for c in self.cols_values if c.startswith('z'))
        assert cols_start_x == cols_start_y == cols_start_z, "Inconsistent force column names."
        return sorted(cols_start_x)


# # //////////////////////////////////////////////////////////////////////////////
