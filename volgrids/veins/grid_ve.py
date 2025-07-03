import numpy as np
import pandas as pd

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class GridVolumetricEnergy(vg.Grid):
    RADIUS_FIX = 0.4
    WIDTH_CYLINDERS = 0.25
    HEIGHT_DISKS = 0.25

    # --------------------------------------------------------------------------
    def __init__(self, ms: vg.MolecularSystem, df: pd.DataFrame, kind: str):
        super().__init__(ms)
        self.df = df[df["kind"] == kind].copy()
        self.kind = kind


    # --------------------------------------------------------------------------
    def populate_grid(self):
        for _,row in self.df.iterrows():
            if   row["energy"]  == 0: continue
            if   row["npoints"] == 2: self._process_2p_interaction(row)
            elif row["npoints"] == 3: self._process_3p_interaction(row)
            elif row["npoints"] == 4: self._process_4p_interaction(row)


    # --------------------------------------------------------------------------
    def _process_2p_interaction(self, row):
        a,b = self._get_positions(row)

        pos = (a + b) / 2  # kernel is placed at the center of the two particles
        direction = b - a  # get the direction vector between two particles

        ##### perform kernel operations
        radius = np.linalg.norm(direction) * self.RADIUS_FIX
        kernel = vg.KernelCylinder(
            radius = radius, vdirection = direction, width = self.WIDTH_CYLINDERS,
            deltas = self.ms.deltas, dtype = np.float32
        )
        self._apply_kernel(kernel, pos, row["energy"])


    # --------------------------------------------------------------------------
    def _process_3p_interaction(self, row):
        a, b, c = self._get_positions(row)
        u = vg.Math.normalize(a - b)
        v = vg.Math.normalize(c - b)
        p = vg.Math.normalize(b - a)
        q = vg.Math.normalize(c - a)

        pos = b                                     # kernel is placed at the vertix B of the triangle ABC
        direction = (u + v) / np.linalg.norm(u + v) # get the direction vector that bisects the angle between the two sides AB and BC
        angle = np.arccos(np.dot(u, v))             # get the angle between the two sides AB and BC
        normal = np.cross(p, q)                     # get the normal vector perpendicular to the plane ABC

        ##### perform kernel operations
        kernel = vg.KernelDiskConecut(
            radius = 2, vnormal = normal, height = self.HEIGHT_DISKS,
            vdirection = direction, max_angle = angle,
            deltas = self.ms.deltas, dtype = np.float32
        )
        self._apply_kernel(kernel, pos, row["energy"])


    # --------------------------------------------------------------------------
    def _process_4p_interaction(self, row):
        a,b,c,d = self._get_positions(row)

        pos = (a + b + c + d) / 4 # kernel is placed at the center of the four particles
        direction = a - c         # get the direction vector between the first  two non-adjacent particles
        normal = b - d            # get the direction vector between the second two non-adjacent particles

        ##### perform kernel operations
        radius0 = np.linalg.norm(direction) * self.RADIUS_FIX
        radius1 = np.linalg.norm(normal) * self.RADIUS_FIX
        kernel0 = vg.KernelCylinder(
            radius = radius0, vdirection = direction, width = self.WIDTH_CYLINDERS,
            deltas = self.ms.deltas, dtype = np.float32
        )
        kernel1 = vg.KernelCylinder(
            radius = radius1, vdirection = normal, width = self.WIDTH_CYLINDERS,
            deltas = self.ms.deltas, dtype = np.float32
        )
        self._apply_kernel(kernel0, pos, row["energy"])
        self._apply_kernel(kernel1, pos, row["energy"])


    # ------------------------------------------------------------------------------
    def _get_positions(self, row: pd.Series):
        ### idxs are expected to be 0-based
        def _split_idx_group(str_idxs: str) -> list[int]:
            return [int(idx) for idx in str_idxs.split('-')]

        if row["idxs_are_residues"]:
            return (
                self.ms.system.residues[idx].atoms.center_of_geometry()
                for idx in _split_idx_group(row["idxs"])
            )

        return self.ms.system.atoms[_split_idx_group(row["idxs"])].positions


    # --------------------------------------------------------------------------
    def _apply_kernel(self, kernel: vg.Kernel, position, energy):
        operation = "max" if energy > 0 else "min"
        kernel.link_to_grid(self.grid, self.ms.minCoords)
        kernel.stamp(position, multiplication_factor = energy, operation = operation)


# //////////////////////////////////////////////////////////////////////////////
