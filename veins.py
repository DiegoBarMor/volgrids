import warnings
import volgrids.veins as ve

if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    ve.AppVeins.from_cli().run()
