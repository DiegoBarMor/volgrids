import math
from dataclasses import dataclass, field
from functools import cached_property, total_ordering
from typing import Dict, List, Optional, Tuple

import numpy
import numpy.typing

from volgrids._vendors.rnapolis.common import (
    Residue,
    ResidueAuth,
    ResidueLabel,
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


def torsion_angle(a1: Atom, a2: Atom, a3: Atom, a4: Atom) -> float:
    """Calculates the torsion angle between four atoms."""
    v1 = a2.coordinates - a1.coordinates
    v2 = a3.coordinates - a2.coordinates
    v3 = a4.coordinates - a3.coordinates

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
