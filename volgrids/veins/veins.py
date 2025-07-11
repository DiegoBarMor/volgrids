import pandas as pd

import volgrids as vg
import volgrids.veins as ve

# ------------------------------------------------------------------------------
def _assert_df(df: pd.DataFrame, *cols_metadata):
    if not set(df.columns).issuperset(cols_metadata):
        raise ValueError(
            f"CSV file '{ve.PATH_ENERGIES_CSV}' must contain the columns: " +\
            ", ".join(map(lambda x: f"'{x}'", cols_metadata)) + " " +\
            f"Found columns: {df.columns}"
        )


# //////////////////////////////////////////////////////////////////////////////
class VeinsApp:
    def __init__(self, *args, **kwargs):
        handler = ve.ParamHandlerVeins(*args, **kwargs)
        handler.assign_globals()

        self.ms = vg.MolSystem(ve.PATH_STRUCTURE, ve.PATH_TRAJECTORY)
        self.df = pd.read_csv(ve.PATH_ENERGIES_CSV).dropna(how = "any")
        self.cols_frames: list = None

        if self.ms.do_traj:
            _assert_df(self.df, "kind", "npoints", "idxs", "idxs_are_residues")
            self.cols_frames = sorted(filter(lambda x: x.startswith("frame"), self.df.columns))
            if not self.cols_frames:
                raise ValueError(
                    f"CSV file '{ve.PATH_ENERGIES_CSV}' must contain at least one column starting with 'frame' "
                    "when running in trajectory mode."
                )
            mat = self.df[self.cols_frames].to_numpy()
            mat[mat < ve.ENERGY_CUTOFF] = 0.0
            self.df.loc[:, self.cols_frames] = mat

        else:
            _assert_df(self.df, "kind", "npoints", "idxs", "idxs_are_residues", "energy")
            self.df = self.df[self.df["energy"].abs() > ve.ENERGY_CUTOFF]

        self.timer = vg.Timer(
            f">>> Now processing '{self.ms.molname}' ({vg.USER_MODE})"
        )


    # --------------------------------------------------------------------------
    @classmethod
    def from_cli(cls):
        params_pos, params_kwd = ve.ParamHandlerVeins.parse_cli_args()
        return cls(*params_pos, **params_kwd)


    # --------------------------------------------------------------------------
    def run(self):
        self.timer.start()

        if self.ms.do_traj: # TRAJECTORY MODE
            for _ in self.ms.system.trajectory:
                current_col = self.cols_frames[self.ms.frame]
                self.df["energy"] = self.df[current_col]
                self.ms.frame += 1

                timer_frame = vg.Timer(f"...>>> Frame {self.ms.frame}/{len(self.ms.system.trajectory)}")
                timer_frame.start()
                self._process_grids()
                timer_frame.end()

        else: # SINGLE PDB MODE
            self._process_grids()

        self.timer.end()


    # --------------------------------------------------------------------------
    def _process_grids(self):
        for kind in self.df["kind"].unique():
            grid = ve.GridVolumetricEnergy(self.ms, self.df, kind)
            grid.populate_grid()
            grid.save_data(ve.FOLDER_OUT, grid.kind)


# # //////////////////////////////////////////////////////////////////////////////
