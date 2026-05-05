import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class AppOccupancy(sm.AppSmiffer):
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main, str_mode = "OGs")
        app_main.load_configs(su)


    # --------------------------------------------------------------------------
    def _process_grids(self):
        if sm.DO_SMIF_HYDROPHOBIC:
            su.OgHydrophobic(self.ms).populate_grid(self.grid_smif)
            sm.Smif.save_data(self.grid_smif, self.ms, self.folder_out, "og_hydrophobic")

        if sm.DO_SMIF_HBA:
            su.OgHBAccepts(self.ms).populate_grid(self.grid_smif)
            sm.Smif.save_data(self.grid_smif, self.ms, self.folder_out, "og_hbacceptors")

        if sm.DO_SMIF_HBD:
            su.OgHBDonors(self.ms).populate_grid(self.grid_smif)
            sm.Smif.save_data(self.grid_smif, self.ms, self.folder_out, "og_hbdonors")

        if sm.DO_SMIF_STACKING:
            su.OgStacking(self.ms).populate_grid(self.grid_smif)
            sm.Smif.save_data(self.grid_smif, self.ms, self.folder_out, "og_stacking")


# //////////////////////////////////////////////////////////////////////////////
