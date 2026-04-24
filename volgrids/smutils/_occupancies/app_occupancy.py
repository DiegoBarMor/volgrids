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
            self._save_og(
                self._calc_smif(su.OgHydrophobic), title = "og_hydrophobic"
            )

        if sm.DO_SMIF_HBA:
            self._save_og(
                self._calc_smif(su.OgHBAccepts), title = "og_hbacceptors"
            )

        if sm.DO_SMIF_HBD:
            self._save_og(
                self._calc_smif(su.OgHBDonors), title = "og_hbdonors"
            )

        if sm.DO_SMIF_STACKING:
            self._save_og(
                self._calc_smif(su.OgStacking), title = "og_stacking"
            )


    # --------------------------------------------------------------------------
    def _save_og(self, smif: sm.Smif,title: str) -> None:
        smif.save_data_smif(smif.grid, self.ms, self.folder_out, title)


# //////////////////////////////////////////////////////////////////////////////
