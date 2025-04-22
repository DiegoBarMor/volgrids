import numpy as np
import volpot as vp

# //////////////////////////////////////////////////////////////////////////////
class KernelSphere(vp.Kernel):
    """For generating simple boolean spheres (e.g. for masks)"""
    def __init__(self, radius, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        self.kernel[self.dist < radius] = 1


# //////////////////////////////////////////////////////////////////////////////
class KernelCylinder(vp.Kernel):
    """For generating boolean cylinders"""
    def __init__(self, radius, vdirection, height, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        projection = vp.get_projection(self.shifted_coords, vdirection)
        norm_coords = vp.get_norm(self.shifted_coords)
        h = np.sqrt(norm_coords**2 - projection**2 )
        self.kernel[h < height] = 1


# //////////////////////////////////////////////////////////////////////////////
class KernelDisk(KernelSphere):
    """For generating boolean disks"""
    def __init__(self, radius, vnormal, height, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        projection = vp.get_projection(self.shifted_coords, vnormal)
        projection = np.abs(projection)
        self.kernel[projection >= height] = 0


# //////////////////////////////////////////////////////////////////////////////
