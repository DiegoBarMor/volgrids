import pandas as pd

import volgrids as vg
import volgrids.veins as ve

# //////////////////////////////////////////////////////////////////////////////
class VeinsApp:
    def __init__(self):
        ve.VeinsArgsParser()

        df = pd.read_csv(ve.PATH_ENERGIES_CSV).dropna(how = "any")

        if not set(df.columns).issuperset({"kind", "npoints", "idxs", "idxs_are_residues", "energy"}):
            raise ValueError(
                f"CSV file '{ve.PATH_ENERGIES_CSV}' must contain the columns: "
                "'kind', 'npoints', 'idxs', 'idxs_are_residues', 'energy'. "
                f"Found columns: {df.columns}"
            )

        self.df = df[df["energy"].abs() > ve.ENERGY_CUTOFF].copy()
        self.ms = vg.MolecularSystem(ve.PATH_STRUCTURE)
        self.timer = vg.Timer(
            f">>> Now processing '{self.ms.molname}' ({vg.USER_MODE})"
        )


    # --------------------------------------------------------------------------
    def run(self):
        self.timer.start()

        if self.ms.do_traj: # TRAJECTORY MODE
            raise NotImplementedError()
            # for _ in self.ms.system.trajectory:
            #     self.ms.frame += 1
            #     timer_frame = vg.Timer(f"...>>> Frame {self.ms.frame}/{len(self.ms.system.trajectory)}")
            #     timer_frame.start()
            #     self._process_grids()
            #     timer_frame.end()

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
