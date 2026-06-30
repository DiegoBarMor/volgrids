import volgrids as vg
import volgrids.smiffer as smf
import volgrids.smutils as sut

# //////////////////////////////////////////////////////////////////////////////
class AppOccupancy(smf.AppSmiffer):
    EXTENSION = ".og" # occupancy grid

    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main, str_mode = "OGs")
        app_main.load_configs()


    # --------------------------------------------------------------------------
    def _process_grids(self):
        self.grid_smif = vg.Grid(self.ms.box, init_grid = not vg.CFG.smif_apbs)

        if vg.CFG.smif_hphob:
            sut.OgHydrophobic(self.ms).populate_grid(self.grid_smif)
            path_out, key_out = self.paths_out["hphob"], self.keys_out["hphob"]
            smf.Smif.save_data(self.grid_smif, self.ms, path_out, key_out)

        if vg.CFG.smif_hba:
            sut.OgHBAccepts(self.ms).populate_grid(self.grid_smif)
            path_out, key_out = self.paths_out["hba"], self.keys_out["hba"]
            smf.Smif.save_data(self.grid_smif, self.ms, path_out, key_out)

        if vg.CFG.smif_hbd:
            sut.OgHBDonors(self.ms).populate_grid(self.grid_smif)
            path_out, key_out = self.paths_out["hbd"], self.keys_out["hbd"]
            smf.Smif.save_data(self.grid_smif, self.ms, path_out, key_out)

        if vg.CFG.smif_stk:
            sut.OgStacking(self.ms).populate_grid(self.grid_smif)
            path_out, key_out = self.paths_out["stk"], self.keys_out["stk"]
            smf.Smif.save_data(self.grid_smif, self.ms, path_out, key_out)


# //////////////////////////////////////////////////////////////////////////////
