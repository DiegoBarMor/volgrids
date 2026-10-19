import math
from dataclasses import dataclass, field
from functools import cached_property, total_ordering
from typing import Dict, List, Optional, Tuple

import numpy
import numpy.typing

from volgrids._vendors.rnapolis.common import (
    BaseInteractions,
    Residue,
    ResidueAuth,
    ResidueLabel,
    Structure2D,
)

BASE_ATOMS = {
    "A": ["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"],
    "G": ["N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"],
    "C": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
    "U": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"],
    "T": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "C7"],
}

BASE_DONORS = {
    "A": ["C2", "N6", "C8", "O2'"],
    "G": ["N1", "N2", "C8", "O2'"],
    "C": ["N4", "C5", "C6", "O2'"],
    "U": ["N3", "C5", "C6", "O2'"],
    "T": ["N3", "C6", "C7"],
}

BASE_ACCEPTORS = {
    "A": ["N1", "N3", "N7"],
    "G": ["N3", "O6", "N7"],
    "C": ["O2", "N3"],
    "U": ["O2", "O4"],
    "T": ["O2", "O4"],
}

PHOSPHATE_ACCEPTORS = ["OP1", "OP2", "O5'", "O3'"]

RIBOSE_ACCEPTORS = ["O4'", "O2'"]

BASE_EDGES = {
    "A": {
        "N1": "W",
        "C2": "WS",
        "N3": "S",
        "N6": "WH",
        "N7": "H",
        "C8": "H",
        "O2'": "S",
    },
    "G": {
        "N1": "W",
        "N2": "WS",
        "N3": "S",
        "O6": "WH",
        "N7": "H",
        "C8": "H",
        "O2'": "S",
    },
    "C": {
        "O2": "WS",
        "N3": "W",
        "N4": "WH",
        "C5": "H",
        "C6": "H",
        "O2'": "S",
    },
    "U": {
        "O2": "WS",
        "N3": "W",
        "O4": "WH",
        "C5": "H",
        "C6": "H",
        "O2'": "S",
    },
    "T": {
        "O2": "WS",
        "N3": "W",
        "O4": "WH",
        "C6": "H",
        "C7": "H",
    },
}


@dataclass(frozen=True, order=True)
class Atom:
    """Representation of a single atom with 3D coordinates.

    Attributes:
        entity_id: Optional entity identifier from mmCIF.
        label: Label-style residue identifier (label_asym_id, label_seq_id, etc.).
        auth: Author-style residue identifier (auth_asym_id, auth_seq_id, etc.).
        model: Model number in the structure (for multi-model files).
        name: Atom name (e.g. "C1'", "N1").
        x, y, z: Cartesian coordinates in Ångströms.
        occupancy: Optional occupancy value (None if not present).
    """

    entity_id: Optional[str]
    label: Optional[ResidueLabel]
    auth: Optional[ResidueAuth]
    model: int
    name: str
    x: float
    y: float
    z: float
    occupancy: Optional[float]

    @cached_property
    def coordinates(self) -> numpy.typing.NDArray[numpy.floating]:
        """Return atom coordinates as a NumPy vector [x, y, z]."""
        return numpy.array([self.x, self.y, self.z])


@dataclass(frozen=True)
@total_ordering
class Residue3D(Residue):
    """3D nucleotide residue with helper methods for geometry and classification.

    This class extends the 2D Residue with:

    - model id and one-letter code,
    - full list of Atom objects,
    - detection of glycosidic torsion (chi) and its class (syn/anti),
    - classification as nucleotide based on phosphate/sugar/base features,
    - utilities for getting base normals, connectivity, and representative atoms.
    """

    model: int
    one_letter_name: str
    atoms: Tuple[Atom, ...]
    standard_residue_name: str = ""

    def __post_init__(self):
        """Default standard_residue_name to self.name when not provided."""
        if not self.standard_residue_name:
            object.__setattr__(self, "standard_residue_name", self.name or "")

    # Dict representing expected name of atom involved in glycosidic bond
    outermost_atoms = {"A": "N9", "G": "N9", "C": "N1", "U": "N1", "T": "N1"}
    # Dist representing expected name of atom closest to the tetrad center
    innermost_atoms = {"A": "N6", "G": "O6", "C": "N4", "U": "O4", "T": "O4"}
    # Heavy atoms in phosphate and ribose
    phosphate_atoms = {"P", "OP1", "OP2", "O3'", "O5'"}
    sugar_atoms = {"C1'", "C2'", "C3'", "C4'", "C5'", "O4'"}
    # Heavy atoms for each main nucleobase
    nucleobase_heavy_atoms = {
        "A": set(["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"]),
        "G": set(["N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"]),
        "C": set(["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]),
        "U": set(["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]),
        "T": set(["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C5M", "C6"]),
    }

    def __lt__(self, other):
        """Order residues by (model, chain, residue number, insertion code)."""
        return (self.model, self.chain or "", self.number, self.icode or " ") < (
            other.model,
            other.chain or "",
            other.number,
            other.icode or " ",
        )

    def __hash__(self):
        """Hash residue by model and residue identifiers."""
        return hash((self.model, self.label, self.auth))

    @cached_property
    def is_nucleotide(self) -> bool:
        """Heuristically determine whether this residue behaves like a nucleotide.

        The decision is based on:

        - presence of phosphate atoms,
        - presence of sugar atoms,
        - match to expected nucleobase heavy atoms,
        - connectivity scores between sugar and base/phosphate.

        Returns True if the combined score exceeds 0.5.
        """
        scores = {"phosphate": 0.0, "sugar": 0.0, "base": 0.0, "connections": 0.0}
        weights = {"phosphate": 0.25, "sugar": 0.25, "base": 0.25, "connections": 0.25}

        residue_atoms = {atom.name for atom in self.atoms}

        phosphate_match = len(residue_atoms.intersection(self.phosphate_atoms))
        scores["phosphate"] = phosphate_match / len(self.phosphate_atoms)

        sugar_match = len(residue_atoms.intersection(self.sugar_atoms))
        scores["sugar"] = sugar_match / len(self.sugar_atoms)

        nucleobase_atoms = {
            key: self.nucleobase_heavy_atoms[key] for key in self.nucleobase_heavy_atoms
        }
        matches = {
            key: len(residue_atoms.intersection(nucleobase_atoms[key]))
            / len(nucleobase_atoms[key])
            for key in nucleobase_atoms
        }
        best_match = max(matches.items(), key=lambda x: x[1])
        scores["base"] = best_match[1]

        connection_score = 0.0
        distance_threshold = 2.0

        if "P" in residue_atoms and "O5'" in residue_atoms:
            p_atom = next(atom for atom in self.atoms if atom.name == "P")
            o5_atom = next(atom for atom in self.atoms if atom.name == "O5'")
            if (
                numpy.linalg.norm(p_atom.coordinates - o5_atom.coordinates)
                <= distance_threshold
            ):
                connection_score += 0.5
        if "C1'" in residue_atoms:
            c1_atom = next(atom for atom in self.atoms if atom.name == "C1'")
            for base_connection in ["N9", "N1"]:
                if base_connection in residue_atoms:
                    base_atom = next(
                        atom for atom in self.atoms if atom.name == base_connection
                    )
                    if (
                        numpy.linalg.norm(c1_atom.coordinates - base_atom.coordinates)
                        <= distance_threshold
                    ):
                        connection_score += 0.5
                        break

        scores["connections"] = connection_score

        probability = sum(
            scores[component] * weights[component] for component in scores.keys()
        )
        return probability > 0.5

    @cached_property
    def base_normal_vector(self) -> Optional[numpy.typing.NDArray[numpy.floating]]:
        """Return a unit vector normal to the base plane, or None if undefined."""
        if self.one_letter_name in "AG":
            n9 = self.find_atom("N9")
            n7 = self.find_atom("N7")
            n3 = self.find_atom("N3")
            if n9 is None or n7 is None or n3 is None:
                return None
            v1 = n7.coordinates - n9.coordinates
            v2 = n3.coordinates - n9.coordinates
        else:
            n1 = self.find_atom("N1")
            c4 = self.find_atom("C4")
            o2 = self.find_atom("O2")
            if n1 is None or c4 is None or o2 is None:
                return None
            v1 = c4.coordinates - n1.coordinates
            v2 = o2.coordinates - n1.coordinates
        normal: numpy.typing.NDArray[numpy.floating] = numpy.cross(v1, v2)
        return normal / numpy.linalg.norm(normal)

    def find_atom(self, atom_name: str) -> Optional[Atom]:
        """Find an atom with the given name in this residue, or return None."""
        for atom in self.atoms:
            if atom.name == atom_name:
                return atom
        return None


@dataclass
class Structure3D:
    """Container for a 3D structure composed of Residue3D objects.

    Attributes:
        residues: List of Residue3D objects in the structure.
        sequence_by_entity: Mapping from entity ID to the full one-letter
            sequence from the file header, if present.
    """

    residues: List[Residue3D]
    sequence_by_entity: Dict[str, str] = field(default_factory=dict)

    def extract_secondary_structure(
        self,
        base_interactions: BaseInteractions,
    ) -> Structure2D:
        """
        Create a secondary structure representation.

        Args:
            base_interactions: Interactions

        Returns:
            The Structure2D object.
        """
        structure2d = Structure2D(
            base_interactions.base_pairs,
            base_interactions.stackings,
            base_interactions.base_ribose_interactions,
            base_interactions.base_phosphate_interactions,
        )
        return structure2d


def torsion_angle(a1: Atom, a2: Atom, a3: Atom, a4: Atom) -> float:
    """Calculates the torsion angle between four atoms."""
    return calculate_torsion_angle_coords(
        a1.coordinates, a2.coordinates, a3.coordinates, a4.coordinates
    )


def calculate_torsion_angle_coords(
    p1: numpy.typing.NDArray[numpy.floating],
    p2: numpy.typing.NDArray[numpy.floating],
    p3: numpy.typing.NDArray[numpy.floating],
    p4: numpy.typing.NDArray[numpy.floating],
) -> float:
    """Calculates the torsion angle between four points defined by their coordinates."""
    v1 = p2 - p1
    v2 = p3 - p2
    v3 = p4 - p3

    # Normalize vectors to avoid issues with very short vectors
    v1_norm = v1 / numpy.linalg.norm(v1) if numpy.linalg.norm(v1) > 1e-6 else v1
    v2_norm = v2 / numpy.linalg.norm(v2) if numpy.linalg.norm(v2) > 1e-6 else v2
    v3_norm = v3 / numpy.linalg.norm(v3) if numpy.linalg.norm(v3) > 1e-6 else v3

    t1 = numpy.cross(v1_norm, v2_norm)
    t2 = numpy.cross(v2_norm, v3_norm)
    t3 = v1_norm * numpy.linalg.norm(v2_norm)

    # Ensure t1 and t2 are not zero vectors before calculating dot products
    if numpy.linalg.norm(t1) < 1e-6 or numpy.linalg.norm(t2) < 1e-6:
        return 0.0  # Or handle as undefined/error

    dot_t1_t2 = numpy.dot(t1, t2)
    dot_t2_t3 = numpy.dot(t2, t3)

    # Clamp dot product arguments for acos/atan2 to avoid domain errors
    dot_t1_t2 = numpy.clip(dot_t1_t2, -1.0, 1.0)

    angle = math.atan2(dot_t2_t3, dot_t1_t2)
    return angle if not math.isnan(angle) else 0.0
