from dataclasses import dataclass

# //////////////////////////////////////////////////////////////////////////////
@dataclass
class SphereInfo:
    x: float
    y: float
    z: float
    radius: float

    # --------------------------------------------------------------------------
    def get_pos(self) -> tuple[float, float, float]:
        return self.x, self.y, self.z

    # --------------------------------------------------------------------------
    def get_str_query(self, extra_dist: float = 0.0) -> str:
        return f"{self.x} {self.y} {self.z} {self.radius + extra_dist}"


# //////////////////////////////////////////////////////////////////////////////
