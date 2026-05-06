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
        """
        The `only_pdb2pqr` argument allows smiffer to split the APBS calculation in two halves
        (first PQR generation, then APBS itself). The temporary PQR and IN files generated in the
        first half are stored in memory (`vg.TMP_APBS_CONTENT_IN`, `vg.TMP_APBS_CONTENT_PQR`) for the
        second half to use. This is controlled by `_run_pdb2pqr_if_needed`. Skipping the first half
        (by either passing `only_pdb2pqr=False` or by calling `run_subprocess_apbs` directly) will
        implicitly run the PQR generation first (this is handled internally inside `apbs.sh`)
        """
        self.atoms = atoms
        self.name_pdb = name_pdb
        self.only_pdb2pqr = only_pdb2pqr
        self._tmpdir: tempfile.TemporaryDirectory = None


    # --------------------------------------------------------------------------
    def __enter__(self) -> Path:
        """Note that the path return depends on the `only_pdb2pqr` parameter"""
        self._tmpdir = tempfile.TemporaryDirectory()

        path_tmpdir   = Path(self._tmpdir.name)
        path_tmp_pdb  = path_tmpdir / self.name_pdb
        path_tmp_pqr  = path_tmpdir / f"{self.name_pdb}.pqr"
        path_tmp_apbs = path_tmpdir / f"{self.name_pdb}.dx"

        ##### First half: PQR and IN generation
        self._run_pdb2pqr_if_needed()

        if self.only_pdb2pqr:
            return path_tmp_pqr

        ##### Second half: APBS calculation
        self._assert_success(
            self.run_subprocess_apbs([str(path_tmp_pdb), "--keep-pqr"])
        )

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
    def _run_pdb2pqr_if_needed(self):
        path_tmpdir   = Path(self._tmpdir.name)
        path_tmp_pdb  = path_tmpdir / self.name_pdb
        path_tmp_in   = path_tmpdir / f"{self.name_pdb}.in"
        path_tmp_pqr  = path_tmpdir / f"{self.name_pdb}.pqr"

        ##### APBS is being run after previously running the PQR generation (loaded into memory)
        if vg.TMP_APBS_CONTENT_IN and vg.TMP_APBS_CONTENT_PQR:
            path_tmp_in.write_text(vg.TMP_APBS_CONTENT_IN)
            path_tmp_pqr.write_text(vg.TMP_APBS_CONTENT_PQR)
            return

        ##### PQR is to be generated and loaded into memory
        self.atoms.write(path_tmp_pdb)
        self._assert_success(
            self.run_subprocess_pdb2pqr([str(path_tmp_pdb)])
        )

        vg.TMP_APBS_CONTENT_IN  = path_tmp_in.read_text()
        vg.TMP_APBS_CONTENT_PQR = path_tmp_pqr.read_text()


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
