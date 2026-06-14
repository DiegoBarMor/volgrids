from dataclasses import dataclass

# //////////////////////////////////////////////////////////////////////////////
@dataclass
class SphereInfo:
    x: float
    y: float
    z: float
    radius: float

    # --------------------------------------------------------------------------
    @classmethod
    def sphere_list(cls, spheres_flat: list[float]) -> list["SphereInfo"]:
        if not spheres_flat: return []
        if len(spheres_flat) % 4: raise ValueError(
            f"Spheres should be provided as a list of floats, with 4 floats per sphere (x,y,z,radius). "+\
            f"Got {len(spheres_flat)} floats, which is not a multiple of 4."
        )
        return [
            cls(x, y, z, radius) for x,y,z,radius in
            zip(*[spheres_flat[i::4] for i in range(4)])
        ]


    # --------------------------------------------------------------------------
    @staticmethod
    def assert_sphere_list(spheres: list["SphereInfo"], nframes: int) -> None:
        if nframes > len(spheres):
            raise ValueError(
                f"Number of spheres provided ({len(spheres)}) should be equal or higher than the number of frames in trajectory ({nframes}). " +\
                "Each sphere should correspond to one frame in the trajectory. "
                "Extra provided spheres will be ignored."
            )


    # --------------------------------------------------------------------------
    def get_pos(self) -> tuple[float, float, float]:
        return self.x, self.y, self.z


    # --------------------------------------------------------------------------
    def get_str_query(self, extra_dist: float = 0.0) -> str:
        return f"{self.x} {self.y} {self.z} {self.radius + extra_dist}"


# //////////////////////////////////////////////////////////////////////////////
