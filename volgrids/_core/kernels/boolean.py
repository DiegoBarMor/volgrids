import numpy as np
import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class KernelSphere(vg.Kernel):
    """For generating simple boolean spheres (e.g. for masks)"""
    def __init__(self, radius, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        self.kernel[self.dist < radius] = 1


# //////////////////////////////////////////////////////////////////////////////
class KernelCylinder(vg.Kernel):
    """For generating boolean cylinders"""
    def __init__(self, radius, vdirection, width, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        w = vg.Math.get_projection_height(self.shifted_coords, vdirection)
        self.kernel[w < width] = 1


# //////////////////////////////////////////////////////////////////////////////
class KernelDisk(KernelSphere):
    """For generating boolean disks"""
    def __init__(self, radius, vnormal, height, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        projection = vg.Math.get_projection(self.shifted_coords, vnormal)
        projection = np.abs(projection)
        self.kernel[projection >= height] = 0


# //////////////////////////////////////////////////////////////////////////////
class KernelDiskConecut(KernelDisk):
    """For generating boolean disks with a cone cut"""
    def __init__(self, radius, vnormal, height, vdirection, max_angle, deltas, dtype):
        super().__init__(radius, vnormal, height, deltas, dtype)
        angle = vg.Math.get_angle(self.shifted_coords, vdirection, in_degrees = False)
        self.kernel[angle > max_angle/2] = 0


# //////////////////////////////////////////////////////////////////////////////
