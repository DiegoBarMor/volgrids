import warnings
import volgrids.smiffer as sm

if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    sm.SmifferApp.from_cli().run()
