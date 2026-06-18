import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSMUtils(vg.AppSubcommand):
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        self.func_operation: callable = None


    # --------------------------------------------------------------------------
    def run(self):
        operation = self.main.subcommands.pop(0)
        if operation == "res_nobp" : return self._run_res_nobp()
        if operation == "res_nostk": return self._run_res_nostk()
        if operation == "chemgen"  : return self._run_chemgen()
        if operation == "sphere"   : return self._run_sphere()
        if operation == "occupancy": return self._run_occupancy()
        if operation == "pwoverlap": return self._run_pwoverlap()
        if operation == "box_dim"  : return self._run_box_dim()
        if operation == "log_apbs" : return self._run_log_apbs()
        raise ValueError(f"Unknown operation: {operation}")


    # --------------------------------------------------------------------------
    def _run_res_nobp(self) -> None:
        path_in = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        print(su.ResiduesNucleic.get_residues_nobp(path_in))


    # --------------------------------------------------------------------------
    def _run_res_nostk(self) -> None:
        path_in = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        print(su.ResiduesNucleic.get_residues_nostk(path_in))


    # --------------------------------------------------------------------------
    def _run_chemgen(self) -> None:
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out",
            assertion = fy.PathAssertion.FILE_OUT, default = path_in.with_suffix(".chem")
        )
        path_out.write_text(
            su.ChemTableLigand(path_in).gen_chemtable()
        )


    # --------------------------------------------------------------------------
    def _run_sphere(self) -> None:
        su.AppSpheres(self.main).run()


    # --------------------------------------------------------------------------
    def _run_occupancy(self) -> None:
        su.AppOccupancy(self.main).run()


    # --------------------------------------------------------------------------
    def _run_pwoverlap(self) -> None:
        su.AppPwOverlap(self.main).run()


    # --------------------------------------------------------------------------
    def _run_box_dim(self) -> None:
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)

        box = sm.MolSystem(path_in).box
        xmin, ymin, zmin = box.min_coords
        xmax, ymax, zmax = box.max_coords
        print(f"{xmin:.3f} {xmax:.3f} {ymin:.3f} {ymax:.3f} {zmin:.3f} {zmax:.3f}")


    # --------------------------------------------------------------------------
    def _run_log_apbs(self) -> None:
        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        grid = vg.Grid.load(path_in)
        sm.SmifAPBS.apply_logabs_transform(grid)
        grid.save(path_out)


# //////////////////////////////////////////////////////////////////////////////
