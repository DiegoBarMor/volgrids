import numpy as np
from dataclasses import dataclass

import volgrids as vg
from volgrids._vendors import molsimple as ms

# //////////////////////////////////////////////////////////////////////////////
@dataclass
class SphereInfo:
    x: float
    y: float
    z: float
    radius: float

    # --------------------------------------------------------------------------
    @classmethod
    def parse_sphere_infos(cls, spheres_flat: list[str]) -> list["SphereInfo"]:
        if not spheres_flat: return []

        spheres = vg.NumberLists(spheres_flat, 4)
        if spheres.error: raise ValueError(
            f"Spheres should be provided as a list of floats, with 4 floats per sphere (x,y,z,radius). {spheres.error}"
        )
        return [cls(x, y, z, radius) for x,y,z,radius in spheres.values]


    # --------------------------------------------------------------------------
    @staticmethod
    def assert_sphere_infos(spheres: list["SphereInfo"], nframes: int) -> None:
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
    def filter_particles(self,
        particles: ms.ParticleGroup, extra_dist: float = 0.0
    ) -> ms.ParticleGroup:
        coords = particles.get_positions_numpy()
        dists = np.linalg.norm(coords - np.array(self.get_pos()), axis = 1)
        mask = dists <= (self.radius + extra_dist)
        return particles.select_mask(mask)


# //////////////////////////////////////////////////////////////////////////////
