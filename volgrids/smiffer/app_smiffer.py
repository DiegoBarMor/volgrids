import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as smf
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSmiffer(vg.AppSubcommand):
    EXTENSION = "" # optional extension for derived classes, should start with dot e.g. ".og"

    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain", str_mode: str = "SMIFs"):
        super().__init__(app_main)
        self.str_mode = str_mode # just for printing

        ##### initialized once in __init__
        self.mm: smf.MoleculeManager
        self.timer: vg.Timer
        self.folder_out: Path
        self.path_traj: Path
        self.nproc: int

        self.paths_out: dict[str, Path]
        self.keys_out: dict[str, str]

        #### set in every call of _process_grids
        self.trimmer: smf.Trimmer
        self.grid_smif: vg.Grid

        #### CLI arguments
        smf.PATH_STRUCT = self.main.get_arg_path(
            "path_in",   assertion = fy.PathAssertion.FILE_IN
        )
        smf.PATH_APBS = self.main.get_arg_path(
            "path_apbs", assertion = fy.PathAssertion.FILE_IN,
            allow_none = True
        )
        self.path_traj = self.main.get_arg_path(
            "path_traj", assertion = fy.PathAssertion.FILE_IN,
            allow_none = True
        )
        smf.PATH_CHEM_LIGAND = self.main.get_arg_path(
            "path_chem", assertion = fy.PathAssertion.FILE_IN,
            allow_none = True
        )
        self.folder_out = self.main.get_arg_path(
            "folder_out", assertion = fy.PathAssertion.DIR_OUT,
            default = smf.PATH_STRUCT.parent
        )
        self.nproc = max(1, self.main.get_arg_int("nproc", default = 1))
        self.do_pack_output = self.main.get_arg_bool("pack")

        self._handle_params_configs()
        self._handle_params_resids()
        self._handle_params_sphere()
        self._handle_params_box()
        self._assert_traj_apbs()

        app_main.load_configs()
        self.init_params()

        ### these must be initialized after the configs are loaded,
        ### so it's appropriate to place them at the end of init
        self.mm = smf.MoleculeManager(smf.PATH_STRUCT, self.path_traj)
        self.timer = vg.Timer(
            f">>> {'Sphere' if self.mm.do_use_sphere else 'Whole'} "+\
            f"{fy.Color.magenta(self.str_mode)} for '{fy.Color.yellow(self.mm.molname)}'"
        )

        self.paths_out, self.keys_out = self._get_paths_keys_out(self.folder_out)
        self.mm.enforce_cmap_output |= self.do_pack_output


    # --------------------------------------------------------------------------
    @staticmethod
    def init_params():
        smf.PARAMS_HPHOB = vg.ParamsGaussianUnivariate(
            mu = vg.CFG.param_hphob_dist_mu, sigma = vg.CFG.param_hbhob_dist_sigma,
        )
        smf.PARAMS_HPHIL = vg.ParamsGaussianUnivariate(
            mu = vg.CFG.param_hphil_dist_mu, sigma = vg.CFG.param_hphil_dist_sigma,
        )
        smf.PARAMS_HBA = vg.ParamsGaussianBivariate(
            mu_0 = vg.CFG.param_hba_angle_mu, mu_1 = vg.CFG.param_hba_dist_mu,
            cov_00 = vg.CFG.param_hba_angle_sigma**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = vg.CFG.param_hba_dist_sigma**2,
        )
        smf.PARAMS_HBD_FREE = vg.ParamsGaussianBivariate(
            mu_0 = vg.CFG.param_hbd_free_angle_mu, mu_1 = vg.CFG.param_hbd_free_dist_mu,
            cov_00 = vg.CFG.param_hbd_free_angle_sigma**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = vg.CFG.param_hbd_free_dist_sigma**2,
        )
        smf.PARAMS_HBD_FIXED = vg.ParamsGaussianBivariate(
            mu_0 = vg.CFG.param_hbd_fixed_angle_mu, mu_1 = vg.CFG.param_hbd_fixed_dist_mu,
            cov_00 = vg.CFG.param_hbd_fixed_angle_sigma**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = vg.CFG.param_hbd_fixed_dist_sigma**2,
        )
        smf.PARAMS_STACK = vg.ParamsGaussianBivariate(
            mu_0 = vg.CFG.param_stk_angle_mu, mu_1 = vg.CFG.param_stk_dist_mu,
            cov_00 = vg.CFG.param_stk_cov00, cov_01 = vg.CFG.param_stk_cov01,
            cov_10 = vg.CFG.param_stk_cov10, cov_11 = vg.CFG.param_stk_cov11,
        )

        ### square root of the DIST contribution to smf.COV_STACKING,
        smf.SIGMA_DIST_STACKING = np.sqrt(vg.CFG.param_stk_cov11)


    # --------------------------------------------------------------------------
    def run(self):
        def _end():
            self.timer.end(text = fy.Color.green("volgrids"), minus = smf.APBS_ELAPSED_TIME)
            vg.Utils.delete_traj_locks(self.path_traj)
            if vg.TMP_APBS_CONTENT_PQR:
                path_pqr = self.folder_out / f"{self.mm.molname}.pqr"
                path_pqr.write_text(vg.TMP_APBS_CONTENT_PQR)


        if (smf.PATH_CHEM_LIGAND is not None) and vg.CFG.smif_apbs:
            vg.CFG.smif_apbs = False
            print(f"\n...--- ligand: {fy.Color.red('skipping APBS')} SMIF calculation.", end = ' ', flush = True)

        for path in self.paths_out.values(): # pre-clear stale CMAP outputs once
            vg.GridIO.remove(path)

        self.timer.start()

        ##### 0) SINGLE PDB MODE
        if not self.mm.do_traj:
            self._process_grids()
            return _end()
        #####

        self._warn_chimerax_incompatible()

        print()

        ### 1.a) TRAJECTORY MODE (multiprocessing)
        if self.nproc > 1:
            smf.TrajMultiprocess(self).run(self.mm.nframes)
            return _end()

        ### 1.b) TRAJECTORY MODE (single process)
        for i in range(self.mm.nframes):
            self.process_frame(i)
        return _end()


    # --------------------------------------------------------------------------
    def process_frame(self, frame_idx: int) -> None:
        """Set per-frame state and run `_process_grids`."""
        self.mm.switch_frame(frame_idx)
        timer = vg.Timer(f"...>>> Frame {self.mm.frame}/{self.mm.nframes}")
        timer.start()
        self._process_grids()
        timer.end()


    # --------------------------------------------------------------------------
    def _process_grids(self):
        def run_with_trim_large():
            """Note that APBS should always be the first SMIF executed (in case it needs to instantiate `self.grid_smif`'s internal array)"""
            if vg.CFG.smif_apbs:
                smf.SmifAPBS(mm).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "large")
                path_out, key_out = self.paths_out["apbs"], self.keys_out["apbs"]
                smf.Smif.save_data(self.grid_smif, mm, path_out, key_out)


        def run_with_trim_mid():
            if vg.CFG.smif_hphob:
                smf.SmifHydrophobic(mm).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "mid")
                path_out, key_out = self.paths_out["hphob"], self.keys_out["hphob"]
                smf.Smif.save_data(self.grid_smif, mm, path_out, key_out)

            if vg.CFG.smif_hba:
                smf.SmifHBAccepts(mm).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "mid")
                path_out, key_out = self.paths_out["hba"], self.keys_out["hba"]
                smf.Smif.save_data(self.grid_smif, mm, path_out, key_out)

            if vg.CFG.smif_hbd:
                smf.SmifHBDonors(mm).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "mid")
                path_out, key_out = self.paths_out["hbd"], self.keys_out["hbd"]
                smf.Smif.save_data(self.grid_smif, mm, path_out, key_out)

            if vg.CFG.smif_stk:
                smf.SmifStacking(mm).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "mid")
                path_out, key_out = self.paths_out["stk"], self.keys_out["stk"]
                smf.Smif.save_data(self.grid_smif, mm, path_out, key_out)

            if vg.CFG.trim_save:
                self.trimmer.run_for_saving("mid")
                mask = self.trimmer.get_mask()
                if mask is None:
                    print(fy.Color.red("WARNING: ") + "Trimming mask was requested to be saved, but it is None. Skipping.")
                else:
                    reverse = vg.Grid.reverse(mask) # save the points that are NOT trimmed
                    path_out, key_out = self.paths_out["trim"], self.keys_out["trim"]
                    smf.Smif.save_data(reverse, mm, path_out, key_out)


        def run_with_trim_small():
            if vg.CFG.smif_hphil:
                smf.SmifHydrophilic(mm).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "small")
                path_out, key_out = self.paths_out["hphil"], self.keys_out["hphil"]
                smf.Smif.save_data(self.grid_smif, mm, path_out, key_out)


        mm = self._current_mol_system()
        self.trimmer = smf.Trimmer(mm)

        self.grid_smif = vg.Grid(mm.box, init_grid = not vg.CFG.smif_apbs)

        if self.trimmer.should_do_trim_large(): run_with_trim_large()
        if self.trimmer.should_do_trim_mid()  : run_with_trim_mid()
        if self.trimmer.should_do_trim_small(): run_with_trim_small()

        if vg.CFG.cav_save and self.trimmer.should_do_cavities():
            self.trimmer.run_for_saving("mid")
            path_out, key_out = self.paths_out["cav"], self.keys_out["cav"]
            smf.Smif.save_data(self.trimmer.cavfinder.grid, mm, path_out, key_out)

        vg.TMP_APBS_CONTENT_PQR = "" # clear temporary PQR so that (optionally) subsequent frames recalculate it

        del self.grid_smif # just in case
        del self.trimmer


    # --------------------------------------------------------------------------
    def _get_paths_keys_out(self, folder_out: Path) -> tuple[dict[str, Path], dict[str, str]]:
        """Return two dictionaries: `{kind: path_out}`, `{kind: key_cmap}` for each SMIF kind that is enabled."""
        def _path_key_out(kind: str) -> tuple[Path, str]:
            if self.mm.do_traj:
                path_out = folder_out / f"{self.mm.molname}.{kind}{self.EXTENSION}.cmap"
                key_cmap = self.mm.molname # frame number is to be appended in every iteration of the traj
                return path_out, key_cmap

            if self.do_pack_output: # --pack flag disregards OUT_FORMAT and uses CMAP
                path_out = folder_out / f"{self.mm.molname}.all{self.EXTENSION}.cmap"
                key_cmap = f"{self.mm.molname}.{kind}"
                return path_out, key_cmap

            fmt = vg.GridFormat.from_str(vg.CFG.out_format)
            path_out = folder_out / f"{self.mm.molname}.{kind}{self.EXTENSION}.{fmt.suffix()}"
            return path_out, kind

        paths = {}; keys = {}
        if vg.CFG.smif_hba:   paths["hba"],   keys["hba"]   = _path_key_out("hba") # old extension: hbacceptors
        if vg.CFG.smif_hbd:   paths["hbd"],   keys["hbd"]   = _path_key_out("hbd") # old extension: hbdonors
        if vg.CFG.smif_stk:   paths["stk"],   keys["stk"]   = _path_key_out("stk") # old extension: stacking
        if vg.CFG.cav_save:   paths["cav"],   keys["cav"]   = _path_key_out("cav") # old extension: cavities
        if vg.CFG.smif_apbs:  paths["apbs"],  keys["apbs"]  = _path_key_out("apbs") # old extension: apbs
        if vg.CFG.trim_save:  paths["trim"],  keys["trim"]  = _path_key_out("trim") # old extension: trimming
        if vg.CFG.smif_hphob: paths["hphob"], keys["hphob"] = _path_key_out("hphob") # old extension: hydrophobic
        if vg.CFG.smif_hphil: paths["hphil"], keys["hphil"] = _path_key_out("hphil") # old extension: hydrophilic
        return paths, keys


    # --------------------------------------------------------------------------
    def _current_mol_system(self) -> "smf.MoleculeManager":
        if not vg.CFG.smif_apbs: return self.mm

        ### When dealing with APBS, a first step of PQR generation should always be performed.
        ### PQR can add hydrogen atoms, which should be reflected in the occupancy trimming
        ### operation or be used by HBond SMIFs.
        smf.SmifAPBS(self.mm).gen_pqr()

        chains = self.mm.get_residue_chains() # size: (nresidues,)
        obj = smf.MoleculeManager.from_pqr_data(vg.TMP_APBS_CONTENT_PQR, self.mm.box, chains)

        smf.MoleculeManager.copy_attrs_except_universe(src = self.mm, dst = obj)
        return obj


    # --------------------------------------------------------------------------
    def _handle_params_configs(self):
        configs = self.main.get_arg_str("configs", is_list = True)
        if not configs: return

        for val in configs:
            if '=' in val:
                vg.STR_CUSTOM_CONFIG += val.replace(' ', '\n') + "\n"
                continue

            val_as_path = Path(val)
            if not val_as_path.exists(): self.main.help_and_exit(1,
                f"The specified config path '{val}' does not exist."
            )
            if val_as_path.is_dir(): self.main.help_and_exit(1,
                f"The specified config path '{val}' is a folder, but a file was expected."
            )
            vg.PATHS_CUSTOM_CONFIG.append(val_as_path)


    # --------------------------------------------------------------------------
    def _handle_params_resids(self):
        def _assert_residue(residue: str) -> int:
            splitted = residue.split('.')
            if len(splitted) != 2: self.main.help_and_exit(1,
                f"Invalid residue '{residue}' provided for --residues option." +\
                "The residues must be in the format \"chain_id.resid\" (e.g. \"A.3 A.4 A.5 B.10\")"
            )
            return residue

        residues = self.main.get_arg_str("residues")
        if not residues: return

        splitted = (x for res in residues for x in res.split())
        smf.CUSTOM_RESIDUES = " ".join(_assert_residue(res) for res in splitted)

        if not smf.CUSTOM_RESIDUES: print(
            fy.Color.red("WARNING: ") + "No valid residues provided for --residues option."
        )


    # --------------------------------------------------------------------------
    def _handle_params_sphere(self):
        spheres_flat = self.main.get_arg_str("sphere", is_list = True)
        if not spheres_flat: return

        try: smf.SPHERES = vg.SphereInfo.parse_sphere_infos(spheres_flat)
        except ValueError as e: self.main.help_and_exit(1, f"{e}")


    # --------------------------------------------------------------------------
    def _handle_params_box(self):
        boxes_flat = self.main.get_arg_str("box", is_list = True)
        if not boxes_flat: return

        if smf.SPHERES: self.main.help_and_exit(1,
            "The --box-csv and --sphere options are mutually exclusive (both define "
            "the grid box). Please provide only one of them."
        )

        try: box_infos = vg.BoxInfo.parse_box_infos(boxes_flat)
        except ValueError as e: self.main.help_and_exit(1, f"{e}")

        smf.BOXES_ENFORCED = [box_info.create_box() for box_info in box_infos]


    # --------------------------------------------------------------------------
    def _warn_chimerax_incompatible(self):
        """Warn when the grid box changes across frames: such a CMAP is a valid HDF5
        file but cannot be opened as a map series by ChimeraX (which requires every
        frame to share the same grid dimensions and origin)."""
        def print_warning(): print(
            fy.Color.red("WARNING:"),
            "the grid box changes across frames, so each frame of the output CMAP",
            "will have a different origin/size. The file is a valid HDF5 CMAP, but",
            fy.Color.red("its visualization is currently not done properly by ChimeraX."),
            "In the case of boxes with different sizes, ChimeraX will not open them as a",
            "map series, instead it willl load all the grids at once. In the case of boxes of same size",
            "but different origins, the grids will be opened as a series but use the box of the first frame.",
        )
        if not self.mm.do_traj: return
        if vg.CFG.box_tight_traj: return print_warning()
        if not smf.BOXES_ENFORCED: return
        if any(
            box != smf.BOXES_ENFORCED[0] for box in smf.BOXES_ENFORCED
        ): print_warning()


    # --------------------------------------------------------------------------
    def _assert_traj_apbs(self):
        if not vg.CFG.smif_apbs: return
        if (self.path_traj is None) or (smf.PATH_APBS is None): return
        self.main.help_and_exit(1,
            f"The APBS output '{smf.PATH_APBS}' was provided. However, "+
            "trajectory mode is enabled, so this file would be ambiguous. "+
            "Please either disable trajectory mode or remove the APBS file input. "
        )


# //////////////////////////////////////////////////////////////////////////////
