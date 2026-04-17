from dataclasses import dataclass

# //////////////////////////////////////////////////////////////////////////////
@dataclass
class SphereInfo:
    x: float
    y: float
    z: float
    radius: float

    # --------------------------------------------------------------------------
    def get_str_query(self, extra_dist: float = 0.0) -> str:
        return f"{self.x} {self.y} {self.z} {self.radius + extra_dist}"


# //////////////////////////////////////////////////////////////////////////////
