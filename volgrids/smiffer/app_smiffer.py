import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSmiffer(vg.AppSubcommand):
    EXTENSION = "" # optional extension for derived classes, should start with dot e.g. ".og"

    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain", str_mode: str = "SMIFs"):
        super().__init__(app_main)
        self.str_mode = str_mode # just for printing

        ##### initialized once in __init__
        self.ms: sm.MolSystem
        self.timer: vg.Timer
        self.folder_out: Path
        self.path_traj: Path
        self.nproc: int

        self.paths_out: dict[str, Path]
        self.keys_out: dict[str, str]

        #### set in every call of _process_grids
        self.trimmer: sm.Trimmer
        self.grid_smif: vg.Grid

        #### CLI arguments
        sm.PATH_STRUCT = self.main.get_arg_path(
            "path_in",   assertion = fy.PathAssertion.FILE_IN
        )
        sm.PATH_APBS = self.main.get_arg_path(
            "path_apbs", assertion = fy.PathAssertion.FILE_IN,
            allow_none = True
        )
        self.path_traj = self.main.get_arg_path(
            "path_traj", assertion = fy.PathAssertion.FILE_IN,
            allow_none = True
        )
        sm.PATH_CHEM_LIGAND = self.main.get_arg_path(
            "path_chem", assertion = fy.PathAssertion.FILE_IN,
            allow_none = True
        )
        self.folder_out = self.main.get_arg_path(
            "folder_out", assertion = fy.PathAssertion.DIR_OUT,
            default = sm.PATH_STRUCT.parent
        )
        self.nproc = max(1, self.main.get_arg_int("nproc", default = 1))
        self.do_pack_output = self.main.get_arg_bool("pack")

        self._handle_params_configs()
        self._handle_params_resids()
        self._handle_params_sphere()
        self._handle_params_box()
        self._assert_traj_apbs()

        app_main.load_configs(vg, sm)
        self.init_params()

        ### these must be initialized after the configs are loaded,
        ### so it's appropriate to place them at the end of init
        self.ms = sm.MolSystem(sm.PATH_STRUCT, self.path_traj)
        self.timer = vg.Timer(
            f">>> {'Sphere' if self.ms.do_ps else 'Whole'} "+\
            f"{fy.Color.magenta(self.str_mode)} for '{fy.Color.yellow(self.ms.molname)}'"
        )

        self.paths_out, self.keys_out = self._get_paths_keys_out(self.folder_out)
        self.ms.enforce_cmap_output |= self.do_pack_output


    # --------------------------------------------------------------------------
    @staticmethod
    def init_params():
        sm.PARAMS_HPHOB = vg.ParamsGaussianUnivariate(
            mu = sm.MU_HYDROPHOBIC, sigma = sm.SIGMA_HYDROPHOBIC,
        )
        sm.PARAMS_HPHIL = vg.ParamsGaussianUnivariate(
            mu = sm.MU_HYDROPHILIC, sigma = sm.SIGMA_HYDROPHILIC,
        )
        sm.PARAMS_HBA = vg.ParamsGaussianBivariate(
            mu_0 = sm.MU_ANGLE_HBA, mu_1 = sm.MU_DIST_HBA,
            cov_00 = sm.SIGMA_ANGLE_HBA**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = sm.SIGMA_DIST_HBA**2,
        )
        sm.PARAMS_HBD_FREE = vg.ParamsGaussianBivariate(
            mu_0 = sm.MU_ANGLE_HBD_FREE, mu_1 = sm.MU_DIST_HBD_FREE,
            cov_00 = sm.SIGMA_ANGLE_HBD_FREE**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = sm.SIGMA_DIST_HBD_FREE**2,
        )
        sm.PARAMS_HBD_FIXED = vg.ParamsGaussianBivariate(
            mu_0 = sm.MU_ANGLE_HBD_FIXED, mu_1 = sm.MU_DIST_HBD_FIXED,
            cov_00 = sm.SIGMA_ANGLE_HBD_FIXED**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = sm.SIGMA_DIST_HBD_FIXED**2,
        )
        sm.PARAMS_STACK = vg.ParamsGaussianBivariate(
            mu_0 = sm.MU_ANGLE_STACKING, mu_1 = sm.MU_DIST_STACKING,
            cov_00 = sm.COV_STACKING_00, cov_01 = sm.COV_STACKING_01,
            cov_10 = sm.COV_STACKING_10, cov_11 = sm.COV_STACKING_11,
        )

        ### square root of the DIST contribution to sm.COV_STACKING,
        sm.SIGMA_DIST_STACKING = np.sqrt(sm.COV_STACKING_11)


    # --------------------------------------------------------------------------
    def run(self):
        def _end():
            self.timer.end(text = fy.Color.green("volgrids"), minus = sm.APBS_ELAPSED_TIME)
            vg.Utils.delete_traj_locks(self.path_traj)
            if vg.TMP_APBS_CONTENT_PQR:
                path_pqr = self.folder_out / f"{self.ms.molname}.pqr"
                path_pqr.write_text(vg.TMP_APBS_CONTENT_PQR)


        if (sm.PATH_CHEM_LIGAND is not None) and sm.DO_SMIF_APBS:
            sm.DO_SMIF_APBS = False
            print(f"\n...--- ligand: {fy.Color.red('skipping APBS')} SMIF calculation.", end = ' ', flush = True)

        for path in self.paths_out.values(): # pre-clear stale CMAP outputs once
            vg.GridIO.remove(path)

        self.timer.start()

        ##### 0) SINGLE PDB MODE
        if not self.ms.do_traj:
            self._process_grids()
            return _end()
        #####

        print()
        n_frames = len(self.ms.system.trajectory)

        ### 1.a) TRAJECTORY MODE (multiprocessing)
        if self.nproc > 1:
            sm.TrajMultiprocess(self).run(n_frames)
            return _end()

        ### 1.b) TRAJECTORY MODE (single process)
        for i in range(n_frames):
            self.process_frame(i)
        return _end()


    # --------------------------------------------------------------------------
    def process_frame(self, frame_idx: int) -> None:
        """Set per-frame state and run `_process_grids`."""
        self.ms.switch_frame(frame_idx)
        n_frames = len(self.ms.system.trajectory)
        timer = vg.Timer(f"...>>> Frame {self.ms.frame}/{n_frames}")
        timer.start()
        self._process_grids()
        timer.end()


    # --------------------------------------------------------------------------
    def _process_grids(self):
        def run_with_trim_large():
            """Note that APBS should always be the first SMIF executed (in case it needs to instantiate `self.grid_smif`'s internal array)"""
            if sm.DO_SMIF_APBS:
                sm.SmifAPBS(ms).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "large")
                path_out, key_out = self.paths_out["apbs"], self.keys_out["apbs"]
                sm.Smif.save_data(self.grid_smif, ms, path_out, key_out)


        def run_with_trim_mid():
            if sm.DO_SMIF_HYDROPHOBIC:
                sm.SmifHydrophobic(ms).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "mid")
                path_out, key_out = self.paths_out["hphob"], self.keys_out["hphob"]
                sm.Smif.save_data(self.grid_smif, ms, path_out, key_out)

            if sm.DO_SMIF_HBA:
                sm.SmifHBAccepts(ms).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "mid")
                path_out, key_out = self.paths_out["hba"], self.keys_out["hba"]
                sm.Smif.save_data(self.grid_smif, ms, path_out, key_out)

            if sm.DO_SMIF_HBD:
                sm.SmifHBDonors(ms).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "mid")
                path_out, key_out = self.paths_out["hbd"], self.keys_out["hbd"]
                sm.Smif.save_data(self.grid_smif, ms, path_out, key_out)

            if sm.DO_SMIF_STACKING:
                sm.SmifStacking(ms).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "mid")
                path_out, key_out = self.paths_out["stk"], self.keys_out["stk"]
                sm.Smif.save_data(self.grid_smif, ms, path_out, key_out)

            if sm.SAVE_TRIMMING_MASK:
                self.trimmer.run_for_saving("mid")
                mask = self.trimmer.get_mask()
                if mask is None:
                    print(fy.Color.red("WARNING: ") + "Trimming mask was requested to be saved, but it is None. Skipping.")
                else:
                    reverse = vg.Grid.reverse(mask) # save the points that are NOT trimmed
                    path_out, key_out = self.paths_out["trim"], self.keys_out["trim"]
                    sm.Smif.save_data(reverse, ms, path_out, key_out)


        def run_with_trim_small():
            if sm.DO_SMIF_HYDROPHILIC:
                sm.SmifHydrophilic(ms).populate_grid(self.grid_smif)
                self.trimmer.trim(self.grid_smif, "small")
                path_out, key_out = self.paths_out["hphil"], self.keys_out["hphil"]
                sm.Smif.save_data(self.grid_smif, ms, path_out, key_out)


        ms = self._current_mol_system()
        self.trimmer = sm.Trimmer(ms)

        self.grid_smif = vg.Grid(ms.box, init_grid = not sm.DO_SMIF_APBS)

        if self.trimmer.should_do_trim_large(): run_with_trim_large()
        if self.trimmer.should_do_trim_mid()  : run_with_trim_mid()
        if self.trimmer.should_do_trim_small(): run_with_trim_small()

        if sm.SAVE_CAVITIES and self.trimmer.should_do_cavities():
            self.trimmer.run_for_saving("mid")
            path_out, key_out = self.paths_out["cav"], self.keys_out["cav"]
            sm.Smif.save_data(self.trimmer.cavfinder.grid, ms, path_out, key_out)

        vg.TMP_APBS_CONTENT_PQR = "" # clear temporary PQR so that (optionally) subsequent frames recalculate it

        del self.grid_smif # just in case
        del self.trimmer


    # --------------------------------------------------------------------------
    def _get_paths_keys_out(self, folder_out: Path) -> tuple[dict[str, Path], dict[str, str]]:
        """Return two dictionaries: `{kind: path_out}`, `{kind: key_cmap}` for each SMIF kind that is enabled."""
        def _path_key_out(kind: str) -> tuple[Path, str]:
            if self.ms.do_traj:
                path_out = folder_out / f"{self.ms.molname}.{kind}{self.EXTENSION}.cmap"
                key_cmap = self.ms.molname # frame number is to be appended in every iteration of the traj
                return path_out, key_cmap

            if self.do_pack_output: # --pack flag disregards GRID_FORMAT_OUTPUT and uses CMAP
                path_out = folder_out / f"{self.ms.molname}.all{self.EXTENSION}.cmap"
                key_cmap = f"{self.ms.molname}.{kind}"
                return path_out, key_cmap

            fmt = vg.GridFormat.from_str(sm.GRID_FORMAT_OUTPUT)
            path_out = folder_out / f"{self.ms.molname}.{kind}{self.EXTENSION}.{fmt.suffix()}"
            return path_out, kind

        paths = {}; keys = {}
        if sm.DO_SMIF_HBA:         paths["hba"],   keys["hba"]   = _path_key_out("hba") # old extension: hbacceptors
        if sm.DO_SMIF_HBD:         paths["hbd"],   keys["hbd"]   = _path_key_out("hbd") # old extension: hbdonors
        if sm.DO_SMIF_STACKING:    paths["stk"],   keys["stk"]   = _path_key_out("stk") # old extension: stacking
        if sm.SAVE_CAVITIES:       paths["cav"],   keys["cav"]   = _path_key_out("cav") # old extension: cavities
        if sm.DO_SMIF_APBS:        paths["apbs"],  keys["apbs"]  = _path_key_out("apbs") # old extension: apbs
        if sm.SAVE_TRIMMING_MASK:  paths["trim"],  keys["trim"]  = _path_key_out("trim") # old extension: trimming
        if sm.DO_SMIF_HYDROPHOBIC: paths["hphob"], keys["hphob"] = _path_key_out("hphob") # old extension: hydrophobic
        if sm.DO_SMIF_HYDROPHILIC: paths["hphil"], keys["hphil"] = _path_key_out("hphil") # old extension: hydrophilic
        return paths, keys


    # --------------------------------------------------------------------------
    def _current_mol_system(self) -> "sm.MolSystem":
        if not sm.DO_SMIF_APBS: return self.ms

        ### When dealing with APBS, a first step of PQR generation should always be performed.
        ### PQR can add hydrogen atoms, which should be reflected in the occupancy trimming
        ### operation or be used by HBond SMIFs.
        sm.SmifAPBS(self.ms).gen_pqr()

        chains = [arr[0] for arr in self.ms.system.residues.chainIDs] # size: (nresidues,)
        obj = sm.MolSystem.from_pqr_data(vg.TMP_APBS_CONTENT_PQR, self.ms.box, chains)

        sm.MolSystem.copy_attributes_except_system(src = self.ms, dst = obj)
        return obj


    # --------------------------------------------------------------------------
    def _handle_params_configs(self):
        configs = self.main.get_arg_str("configs")
        if not configs: return

        if isinstance(configs, bool):
            ### this happens when -c is passed without any value to it
            ### `True` is stored instead of any string value
            available = "\n    " + "\n    ".join(sorted(vg.KNOWN_CONFIGS))
            print(f"Available configuration keys:{available}")
            exit(0)

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
        sm.CUSTOM_RESIDUES = " ".join(_assert_residue(res) for res in splitted)

        if not sm.CUSTOM_RESIDUES: print(
            fy.Color.red("WARNING: ") + "No valid residues provided for --residues option."
        )


    # --------------------------------------------------------------------------
    def _handle_params_sphere(self):
        spheres_flat = self.main.get_arg_float("sphere", is_list = True)
        if not spheres_flat: return

        try: sm.SPHERES = sm.SphereInfo.sphere_list(spheres_flat)
        except ValueError as e: self.main.help_and_exit(1, f"{e}")


    # --------------------------------------------------------------------------
    def _handle_params_box(self):
        box = self.main.get_arg_float("box", is_list = True)
        if not box: return

        x_min, x_max, y_min, y_max, z_min, z_max = box
        min_coords = np.array([x_min, y_min, z_min])
        max_coords = np.array([x_max, y_max, z_max])
        sm.BOX_ENFORCED = vg.Box.from_min_max(min_coords, max_coords)


    # --------------------------------------------------------------------------
    def _assert_traj_apbs(self):
        if not sm.DO_SMIF_APBS: return
        if (self.path_traj is None) or (sm.PATH_APBS is None): return
        self.main.help_and_exit(1,
            f"The APBS output '{sm.PATH_APBS}' was provided. However, "+
            "trajectory mode is enabled, so this file would be ambiguous. "+
            "Please either disable trajectory mode or remove the APBS file input. "
        )


# //////////////////////////////////////////////////////////////////////////////
