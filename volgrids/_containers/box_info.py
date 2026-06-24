from dataclasses import dataclass

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
@dataclass
class BoxInfo:
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    zmin: float
    zmax: float

    # --------------------------------------------------------------------------
    @classmethod
    def from_list(cls, box_list: list[float]) -> "BoxInfo":
        if len(box_list) != 6:
            raise ValueError(
                f"Box should be provided as a list of 6 floats (xmin,xmax,ymin,ymax,zmin,zmax). Provided list has {len(box_list)} elements."
            )
        return cls(*box_list)


    # --------------------------------------------------------------------------
    @classmethod
    def from_box(cls, box: vg.Box) -> "BoxInfo":
        xmin, ymin, zmin = box.min_coords
        xmax, ymax, zmax = box.max_coords
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)


    # --------------------------------------------------------------------------
    def create_box(self) -> "vg.Box":
        return vg.Box.from_min_max(
            [self.xmin, self.ymin, self.zmin],
            [self.xmax, self.ymax, self.zmax],
        )


    # --------------------------------------------------------------------------
    def values(self) -> tuple[float]:
        return self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax


# //////////////////////////////////////////////////////////////////////////////
