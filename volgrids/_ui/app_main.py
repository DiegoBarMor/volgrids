from pathlib import Path

import volgrids as vg
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppMain(fy.App):
    """
    Entry point for the volgrids suite of applications.
    It derives from `freyacli`'s `App` class, which takes care of general CLI handling.
    `AppMain` then dispatches to the specific application (e.g. `AppSmiffer`) based on the first CLI argument, which is expected to be the application's name.
    These specific applications derive from `AppSubcommand`.
    """
    _APP_NAME = "volgrids"
    _VERSION = vg.__version__

    # --------------------------------------------------------------------------
    def __init__(self,  argv: list[str]):
        dir_ui = vg.Utils.resolve_path_package("_ui")
        super().__init__(
            args = argv,
            path_fyr = dir_ui / "fy_rules.fyr",
            path_fyh = dir_ui / "fy_help.fyh",
        )
        self.subcommands = self.get_path_to_root()


    # --------------------------------------------------------------------------
    def run(self):
        app = self.subcommands.pop(0)
        if app == "smiffer": return self._run_smiffer()
        if app == "smutils": return self._run_smutils()
        if app == "vgtools": return self._run_vgtools()
        if app == "apbs"   : return self._run_apbs()


    # --------------------------------------------------------------------------
    def load_configs(self) -> None:
        self._load_config(vg.PATH_DEFAULT_CONFIG, is_file = True)
        for path_config in vg.PATHS_CUSTOM_CONFIG:
            self._load_config(path_config, is_file = True)
        self._load_config(vg.STR_CUSTOM_CONFIG, is_file = False)


    # --------------------------------------------------------------------------
    def _run_smiffer(self) -> None:
        import volgrids.smiffer as sm
        sm.AppSmiffer(self).run()


    # --------------------------------------------------------------------------
    def _run_smutils(self) -> None:
        import volgrids.smutils as su
        su.AppSMUtils(self).run()


    # --------------------------------------------------------------------------
    def _run_vgtools(self) -> None:
        import volgrids.vgtools as vgt
        vgt.AppVGTools(self).run()


    # --------------------------------------------------------------------------
    def _run_apbs(self) -> None:
        ### the parsed flags must be reconstucted.
        ### freyacli is in charge of not letting unexpected flags/arguments through
        cmd = [self.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)]
        if self.get_arg_bool("conv2mrc"): cmd.append("--mrc")
        if self.get_arg_bool("keep_pqr"): cmd.append("--keep-pqr")
        if self.get_arg_bool("verbose" ): cmd.append("--verbose")

        print(f">>> Launching {fy.Color.red('APBS')} subprocess for '{fy.Color.blue(cmd[0])}'...", flush = True)
        apbs = vg.APBSSubprocess.run_subprocess_apbs(cmd)
        print(f"{apbs.stdout}\n{apbs.stderr}".strip(), flush = True)
        exit(apbs.returncode)


    # --------------------------------------------------------------------------
    def _load_config(self, config: Path | str, is_file: bool) -> None:
        if is_file:
            if config is None: return
            ini = vg.ParserIni.from_file(config)
        else:
            if not config.strip(): return
            ini = vg.ParserIni(config)

        try:
            vg.CFG.update_configs_from_ini(ini)
        except ValueError as e:
            self.help_and_exit(1, str(e))


# //////////////////////////////////////////////////////////////////////////////
