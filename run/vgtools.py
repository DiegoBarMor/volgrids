import warnings

### simulate having "volgrids" installed as a package
### this way it's not necessary to install the repo to run this script
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import volgrids.vgtools as vgt

if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    vgt.AppVGTools.from_cli().run()
