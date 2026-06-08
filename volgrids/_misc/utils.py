from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////
class Utils:
    # --------------------------------------------------------------------------
    @staticmethod
    def resolve_path_package(path_to_resolve: str | Path):
        """Resolve the path to the package root directory."""
        root_package = Path(__file__).resolve().parent.parent
        return root_package / path_to_resolve


    # --------------------------------------------------------------------------
    @classmethod
    def assert_vendors(cls):
        """
        Automatic fetching of the vendors first time volgrids is executed,
        in case volgrids isn't downloaded via pip.
        """
        dir_vendors = cls.resolve_path_package("_vendors")
        names_vendors = ("freyacli", "molutils")
        if all((dir_vendors / name).is_dir() for name in names_vendors):
            return

        print(">>> Fetching vendor dependencies for the first time...")
        cls.run_bash(
            "_vendors/fetch_vendors.sh", err_msg = "Failed to fetch vendors"
        )


    # --------------------------------------------------------------------------
    @classmethod
    def run_bash(cls, path_sh: str, *args, err_msg: str = "") -> str:
        import subprocess
        args = (str(arg) for arg in args)
        proc = subprocess.run(
            ["/bin/bash", str(cls.resolve_path_package(path_sh)), *args],
            capture_output = True, text = True
        )
        if proc.returncode == 0: return proc.stdout

        raise RuntimeError('\n'.join((
            f"{err_msg} (code={proc.returncode}):",
            proc.stdout or "<empty_stdout>",
            proc.stderr or "<empty_stderr>",
        )))


    # --------------------------------------------------------------------------
    @staticmethod
    def create_mda_universe_quiet(path_pdb: Path, path_traj: Path|None = None):
        import io
        import logging
        import warnings
        import MDAnalysis as mda
        from contextlib import redirect_stdout, redirect_stderr

        buf = io.StringIO()
        logger = logging.getLogger("MDAnalysis")
        old_level = logger.getEffectiveLevel()
        try:
            logger.setLevel(logging.ERROR)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with redirect_stdout(buf), redirect_stderr(buf):
                    args_traj = [str(path_traj)] if path_traj is not None else []
                    u = mda.Universe(str(path_pdb), *args_traj)
        finally:
            logger.setLevel(old_level)
        return u


    # --------------------------------------------------------------------------
    @staticmethod
    def delete_traj_locks(path_traj: Path|None):
        if path_traj is None: return
        if path_traj.suffix != ".xtc": return

        preffix = str(path_traj.parent / f".{path_traj.stem}.xtc_offsets")
        Path(f"{preffix}.lock").unlink(missing_ok = True)
        Path(f"{preffix}.npz").unlink(missing_ok = True)


# //////////////////////////////////////////////////////////////////////////////
