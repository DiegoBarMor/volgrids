import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class AppOccupancy(sm.AppSmiffer):
    EXTENSION = ".og" # occupancy grid

    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main, str_mode = "OGs")
        app_main.load_configs(su)


    # --------------------------------------------------------------------------
    def _process_grids(self):
        self.grid_smif = vg.Grid(self.ms.box, init_grid = not sm.DO_SMIF_APBS)

        if sm.DO_SMIF_HYDROPHOBIC:
            su.OgHydrophobic(self.ms).populate_grid(self.grid_smif)
            path_out, key_out = self.paths_out["hphob"], self.keys_out["hphob"]
            sm.Smif.save_data(self.grid_smif, self.ms, path_out, key_out)

        if sm.DO_SMIF_HBA:
            su.OgHBAccepts(self.ms).populate_grid(self.grid_smif)
            path_out, key_out = self.paths_out["hba"], self.keys_out["hba"]
            sm.Smif.save_data(self.grid_smif, self.ms, path_out, key_out)

        if sm.DO_SMIF_HBD:
            su.OgHBDonors(self.ms).populate_grid(self.grid_smif)
            path_out, key_out = self.paths_out["hbd"], self.keys_out["hbd"]
            sm.Smif.save_data(self.grid_smif, self.ms, path_out, key_out)

        if sm.DO_SMIF_STACKING:
            su.OgStacking(self.ms).populate_grid(self.grid_smif)
            path_out, key_out = self.paths_out["stk"], self.keys_out["stk"]
            sm.Smif.save_data(self.grid_smif, self.ms, path_out, key_out)


# //////////////////////////////////////////////////////////////////////////////
