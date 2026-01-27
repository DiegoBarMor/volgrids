import tempfile
import subprocess
from pathlib import Path
from MDAnalysis.core.groups import AtomGroup

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class APBSSubprocess:
    _PATH_SCRIPT = vg.resolve_path_package("utils/apbs.sh")

    # --------------------------------------------------------------------------
    def __init__(self, atoms: AtomGroup, name_pdb: str, keep_pqr: bool = False):
        self.atoms = atoms
        self.name_pdb = name_pdb
        self.keep_pqr = keep_pqr
        self._tmpdir = None

    # --------------------------------------------------------------------------
    def __enter__(self) -> Path:
        self._tmpdir = tempfile.TemporaryDirectory()
        path_tmpdir = Path(self._tmpdir.name)

        path_tmp_pdb  = path_tmpdir / self.name_pdb
        path_tmp_apbs = path_tmpdir / f"{self.name_pdb}.dx"

        self.atoms.write(path_tmp_pdb)
        proc = self.run_subprocess([str(path_tmp_pdb), "--pqr"])

        if proc.returncode != 0:
            self._tmpdir.cleanup()
            raise RuntimeError(f"apbs.sh failed (code={proc.returncode}):\n{proc.stderr}")
        if not path_tmp_apbs.exists():
            self._tmpdir.cleanup()
            raise FileNotFoundError(f"Expected APBS output not found: {path_tmp_apbs}")

        if self.keep_pqr:
            vg.PQR_CONTENTS_TEMP = (path_tmpdir / f"{self.name_pdb}.pqr").read_text()

        return path_tmp_apbs

    # --------------------------------------------------------------------------
    def __exit__(self, exc_type, exc, tb):
        if self._tmpdir is not None:
            self._tmpdir.cleanup()
        return False

    # --------------------------------------------------------------------------
    @classmethod
    def run_subprocess(self, args: list[str]) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["/bin/bash", str(self._PATH_SCRIPT)] + args,
            capture_output = True, text = True
        )


# //////////////////////////////////////////////////////////////////////////////
