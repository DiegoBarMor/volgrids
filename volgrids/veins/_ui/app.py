import pandas as pd

import volgrids as vg
import volgrids.veins as ve

DEFAULT_ENERGY_CUTOFF = 1e-3 # [TODO] this should be a config

# ------------------------------------------------------------------------------
def _assert_df(df: pd.DataFrame, *cols_metadata):
    if not set(df.columns).issuperset(cols_metadata):
        raise ValueError(
            f"CSV file '{ve.PATH_ENERGIES_CSV}' must contain the columns: " +\
            ", ".join(map(lambda x: f"'{x}'", cols_metadata)) + " " +\
            f"Found columns: {df.columns}"
        )


# //////////////////////////////////////////////////////////////////////////////
class AppVeins(vg.AppSubcommand):
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        self.func_operation: callable = None


    # --------------------------------------------------------------------------
    def assign_globals(self):
        ve.MODE = self.main.subcommands.pop(0)
        self.func_operation = {
            "energies" : self._run_energies,
            "forces"   : self._run_forces,
        }[ve.MODE]


    # --------------------------------------------------------------------------
    def run(self):
        self.func_operation()


    # --------------------------------------------------------------------------
    def _run_energies(self) -> None:
        self.main.assert_paths(
            keys_file_in = ["path_struct", "path_csv"],
            allow_none = False,
        )
        ve.PATH_STRUCT       = self.main.get_arg_path("path_struct")
        ve.PATH_ENERGIES_CSV = self.main.get_arg_path("path_csv")

        if self.main.get_arg_value("folder_out") is None:
            self.main.set_arg_value("folder_out", ve.PATH_STRUCT.parent)

        self.main.assert_paths(keys_dir_out = ["folder_out"], allow_none = False)
        ve.FOLDER_OUT = self.main.get_arg_path("folder_out")

        ve.PATH_TRAJ = self.main.get_arg_path("path_traj")
        ve.ENERGY_CUTOFF = self.main.get_arg_float("cutoff")
        if ve.ENERGY_CUTOFF is None: ve.ENERGY_CUTOFF = DEFAULT_ENERGY_CUTOFF


        ### [TODO] update this below
        ### init phase
        self.ms = vg.MolSystem(ve.PATH_STRUCT, ve.PATH_TRAJ)
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
            f">>> Now processing '{self.ms.molname}' ({ve.MODE})"
        )

        ### run phase
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
    def _run_forces(self) -> None:
        raise NotImplementedError("VEINS Forces mode is not yet implemented.")


    # --------------------------------------------------------------------------
    def _process_grids(self):
        for kind in self.df["kind"].unique():
            grid = ve.GridVolumetricEnergy(self.ms, self.df, kind)
            grid.populate_grid()
            grid.save_data(ve.FOLDER_OUT, grid.kind)


# # # //////////////////////////////////////////////////////////////////////////////
