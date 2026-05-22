import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppPwOverlap(vg.AppSubcommand):
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)

        self.ms_src: sm.MolSystem
        self.ms_dst: sm.MolSystem
        self.path_out: Path

        path_src = self.main.get_arg_path("path_source", assertion = fy.PathAssertion.FILE_IN)
        path_dst = self.main.get_arg_path("path_target", assertion = fy.PathAssertion.FILE_IN)
        self.path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        app_main.load_configs(vg, sm, su)
        sm.AppSmiffer.init_params()

        self.ms_src = sm.MolSystem(path_src)
        self.ms_dst = sm.MolSystem(path_dst)


    # --------------------------------------------------------------------------
    def run(self):
        smif_src = sm.SmifStacking(self.ms_src)
        smif_dst = sm.SmifStacking(self.ms_dst)

        agroups_dst = list(smif_dst.iter_particles())
        cogs_dst, normals_dst = zip(*(
            sm.SmifStacking.get_cog_normal(ag) for ag in agroups_dst
        ))
        cogs_dst = np.array(cogs_dst).reshape(1, 1, -1, 3)
        normals_dst = np.array(normals_dst).reshape(1, 1, -1, 3)
        arr_dst = np.zeros(cogs_dst.shape[:-1], dtype = vg.FLOAT_DTYPE)

        for agroups_src in smif_src.iter_particles():
            cog_src, normal_src = sm.SmifStacking.get_cog_normal(agroups_src)
            centered_cogs_dst = cogs_dst - cog_src
            dists = vg.Math.get_norm(centered_cogs_dst)

            input_mat = vg.KernelGaussianBivariateAngleDist.get_input_mat(
                centered_cogs_dst, dists, normal_src, is_stacking = True
            )
            arr_dst += vg.Math.bivariate_gaussian(
                input_mat, sm.PARAMS_STACK.mu, sm.PARAMS_STACK.cov_inv
            ) * sm.ENERGY_SCALE


        normals_dst += cogs_dst
        chain_ids = (ag[0].chainID for ag in agroups_dst)
        resids    = (ag[0].resid for ag in agroups_dst)
        residues  = [f"{chain_id}.{resid}" for chain_id, resid in zip(chain_ids, resids)]

        self.path_out.write_text(
            "residue,pwoverlap_stk\n"+
            '\n'.join(
                f"{res},{val:.6f}" for res, val in zip(residues, arr_dst[0][0])
            )
        )


# //////////////////////////////////////////////////////////////////////////////
