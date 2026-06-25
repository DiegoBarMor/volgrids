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
    def parse_box_infos(cls, boxes_flat: list[float]) -> list["BoxInfo"]:
        if not boxes_flat: return []

        boxes = vg.NumberLists(boxes_flat, 6)
        if boxes.error: raise ValueError(
            f"Box should be provided as a list of 6 floats (xmin,xmax,ymin,ymax,zmin,zmax). {boxes.error}"
        )
        return [
            cls(xmin, xmax, ymin, ymax, zmin, zmax)
            for xmin, xmax, ymin, ymax, zmin, zmax in boxes.values
        ]


    # --------------------------------------------------------------------------
    @staticmethod
    def assert_box_infos(boxes: list["BoxInfo"], nframes: int) -> None:
        if nframes > len(boxes):
            raise ValueError(
                f"Number of boxes provided ({len(boxes)}) should be equal or higher than the number of frames in trajectory ({nframes}). " +\
                "Each box should correspond to one frame in the trajectory. "
                "Extra provided boxes will be ignored."
            )


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
