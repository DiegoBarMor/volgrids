import warnings
import volgrids.vgtools as vgt

if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    vgt.AppVGTools.from_cli().run()
