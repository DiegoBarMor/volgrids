import tempfile
import subprocess
from pathlib import Path
from MDAnalysis.core.groups import AtomGroup

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class APBSSubprocess:
    _PATH_SH_APBS    = vg.resolve_path_package("apbs/apbs.sh")
    _PATH_SH_PDB2PQR = vg.resolve_path_package("apbs/pdb2pqr.sh")

    # --------------------------------------------------------------------------
    def __init__(self, atoms: AtomGroup, name_pdb: str, only_pdb2pqr: bool = False):
        self.atoms = atoms
        self.name_pdb = name_pdb
        self.only_pdb2pqr = only_pdb2pqr
        self._tmpdir: tempfile.TemporaryDirectory = None


    # --------------------------------------------------------------------------
    def __enter__(self) -> Path:
        self._tmpdir = tempfile.TemporaryDirectory()
        path_tmpdir = Path(self._tmpdir.name)

        path_tmp_pdb  = path_tmpdir / self.name_pdb
        path_tmp_pqr  = path_tmpdir / f"{self.name_pdb}.pqr"
        path_tmp_apbs = path_tmpdir / f"{self.name_pdb}.dx"

        args = [str(path_tmp_pdb)]

        self.atoms.write(path_tmp_pdb)
        self._assert_success(self.run_subprocess_pdb2pqr(args))

        vg.PQR_CONTENTS_TEMP = path_tmp_pqr.read_text()

        if self.only_pdb2pqr:
            self._safe_cleanup()
            return path_tmp_pqr

        self._assert_success(self.run_subprocess_apbs(args + ["--keep-pqr"]))

        if not path_tmp_apbs.exists():
            self._safe_cleanup()
            raise FileNotFoundError(f"Expected {self._PATH_SH_APBS} output not found: {path_tmp_apbs}")

        return path_tmp_apbs


    # --------------------------------------------------------------------------
    def __exit__(self, exc_type, exc, tb):
        self._safe_cleanup()
        return False


    # --------------------------------------------------------------------------
    @classmethod
    def run_subprocess_pdb2pqr(cls, args: list[str]) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["/bin/bash", str(cls._PATH_SH_PDB2PQR)] + args,
            capture_output = True, text = True
        )


    # --------------------------------------------------------------------------
    @classmethod
    def run_subprocess_apbs(cls, args: list[str]) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["/bin/bash", str(cls._PATH_SH_APBS)] + args,
            capture_output = True, text = True
        )


    # --------------------------------------------------------------------------
    def _assert_success(self, proc: subprocess.CompletedProcess):
        if proc.returncode == 0: return

        self._safe_cleanup()
        raise RuntimeError('\n'.join((
            f"{proc.args[0]} failed (code={proc.returncode}):",
            proc.stderr or "<empty_stderr>",
            proc.stdout or "<empty_stdout>",
        )))

    # --------------------------------------------------------------------------
    def _safe_cleanup(self):
        if self._tmpdir is None: return
        self._tmpdir.cleanup()


# //////////////////////////////////////////////////////////////////////////////
