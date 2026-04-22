from pathlib import Path

import volgrids as vg

try: import freyacli as fy # CLI management
except ImportError: from volgrids._vendors import freyacli as fy

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
        dir_ui = vg.resolve_path_package("_ui")
        super().__init__(
            args = argv,
            path_fyr = dir_ui / "fy_rules.fyr",
            path_fyh = dir_ui / "fy_help.fyh",
        )
        self.subcommands = self.get_path_to_root()

        command = self.subcommands.pop(0)
        func_init = {
            "smiffer": self._init_smiffer,
            "smutils": self._init_smutils,
            "apbs":    self._init_apbs,
            "vgtools": self._init_vgtools,
            "veins":   self._init_veins,
        }[command]
        self.sub_app: "vg.AppSubcommand" = func_init()


    # --------------------------------------------------------------------------
    def run(self):
        self.sub_app.run()


    # --------------------------------------------------------------------------
    def _init_smiffer(self) -> "vg.AppSubcommand":
        import volgrids.smiffer as sm
        app = sm.AppSmiffer(self)
        self._load_all_configs(vg, sm)
        app.init_smif_parameters()
        return app


    # --------------------------------------------------------------------------
    def _init_smutils(self) -> "vg.AppSubcommand":
        import volgrids.smiffer as sm
        import volgrids.smutils as su
        app = su.AppSMUtils(self)
        self._load_all_configs(vg, sm, su)
        return app


    # --------------------------------------------------------------------------
    def _init_apbs(self) -> "vg.AppSubcommand":
        ### the parsed flags must be reconstucted.
        ### freyacli is in charge of not letting unexpected flags/arguments through
        cmd = [self.get_arg_path("path_in")]
        if self.get_arg_bool("conv2mrc"): cmd.append("--mrc")
        if self.get_arg_bool("keep_pqr"): cmd.append("--pqr")
        if self.get_arg_bool("verbose" ): cmd.append("--verbose")

        print(f">>> Launching {fy.Color.red('APBS')} subprocess for '{fy.Color.blue(cmd[0])}'...", flush = True)
        apbs = vg.APBSSubprocess.run_subprocess(cmd)
        print(f"{apbs.stdout}\n{apbs.stderr}".strip(), flush = True)
        exit(apbs.returncode)


    # --------------------------------------------------------------------------
    def _init_vgtools(self) -> "vg.AppSubcommand":
        import volgrids.vgtools as vgt
        app = vgt.AppVGTools(self)
        self._load_all_configs(vg, vgt)
        return app


    # --------------------------------------------------------------------------
    def _init_veins(self) -> "vg.AppSubcommand":
        import volgrids.veins as ve
        app = ve.AppVeins(self)
        self._load_all_configs(vg, ve)
        return app


    # --------------------------------------------------------------------------
    def _load_all_configs(self, *modules: tuple) -> None:
        self._load_config(vg.PATH_DEFAULT_CONFIG, modules)
        for path_config in vg.PATHS_CUSTOM_CONFIG:
            self._load_config(path_config, modules)
        self._load_config(vg.STR_CUSTOM_CONFIG, modules, is_file = False)


    # --------------------------------------------------------------------------
    def _load_config(self,
        config: Path | str, modules: tuple, is_file = True
    ) -> None:
        if is_file:
            if config is None: return
            parser = vg.ParserConfig.from_file(config)
        else:
            if not config.strip(): return
            parser = vg.ParserConfig(config)

        for scope_module in modules:
            parser.apply_config(
                scope_module = scope_module.__dict__,
                this_module_keys = scope_module.__config_keys__,
            )


# //////////////////////////////////////////////////////////////////////////////
