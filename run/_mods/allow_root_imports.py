### simulate having "volgrids" installed as a package
### this way it's not necessary to install the repo to run these scripts
import sys
from pathlib import Path
folder_root = Path(__file__).parent.parent.parent
folder_src = folder_root / "src"
sys.path.insert(0, str(folder_src))
