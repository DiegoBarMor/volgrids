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
            if   row["npoints"] == 2: self._process_2p_interaction(row)
            elif row["npoints"] == 3: self._process_3p_interaction(row)
            elif row["npoints"] == 4: self._process_4p_interaction(row)


    # --------------------------------------------------------------------------
    def _process_2p_interaction(self, row):
        ##### kernel is placed at the center of the two particles
        a,b = self._get_positions(row)
        pos = (a + b) / 2

        ##### get the direction vector between two particles
        a,b = self._get_positions(row)
        direction = b - a

        ##### perform kernel operations
        radius = np.linalg.norm(direction) * self.RADIUS_FIX
        kernel = vg.KernelCylinder(
            radius = radius, vdirection = direction, width = self.WIDTH_CYLINDERS,
            deltas = self.ms.deltas, dtype = np.float32
        )
        self._apply_kernel(kernel, pos, row["energy"])


    # --------------------------------------------------------------------------
    def _process_3p_interaction(self, row):
        ##### kernel is placed at the vertix B of the triangle ABC
        _,b,_ = self._get_positions(row)
        pos = b

        ##### given a triangle ABC, get the direction vector
        ##### that bisects the angle between the two sides AB and BC
        a, b, c = self._get_positions(row)
        u = vg.Math.normalize(a - b)
        v = vg.Math.normalize(c - b)
        direction = (u + v) / np.linalg.norm(u + v)

        ##### get the normal vector perpendicular to the plane ABC
        a, b, c = self._get_positions(row)
        u = vg.Math.normalize(b - a)
        v = vg.Math.normalize(c - a)
        normal = np.cross(u, v)

        ##### given a triangle ABC, get the angle between the two sides AB and BC
        a, b, c = self._get_positions(row)
        u = vg.Math.normalize(a - b)
        v = vg.Math.normalize(c - b)
        angle = np.arccos(np.dot(u, v))

        ##### perform kernel operations
        kernel = vg.KernelDiskConecut(
            radius = 2, vnormal = normal, height = self.HEIGHT_DISKS,
            vdirection = direction, max_angle = angle,
            deltas = self.ms.deltas, dtype = np.float32
        )
        self._apply_kernel(kernel, pos, row["energy"])


    # --------------------------------------------------------------------------
    def _process_4p_interaction(self, row):
        ##### kernel is placed at the center of the four particles
        a,b,c,d = self._get_positions(row)
        pos = (a + b + c + d) / 4

        ##### get the direction vector between the first two non-adjacent particles
        a, _, c, _ = self._get_positions(row)
        direction = a - c

        ##### get the direction vector between the first two non-adjacent particles
        _, b, _, d = self._get_positions(row)
        normal = b - d

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
        def _split_idx_group(str_idxs: str) -> list[int]:
            # idxs are expected to be 0-based
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
