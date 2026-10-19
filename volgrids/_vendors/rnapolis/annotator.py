import math
from collections import Counter
from typing import Dict, List, Optional, Set, Tuple

import numpy
import numpy.typing
from scipy.spatial import KDTree

from volgrids._vendors.rnapolis.common import (
    BaseInteractions,
    BasePair,
    LeontisWesthof,
    Residue,
    Saenger,
    Stacking,
    StackingTopology,
)

from volgrids._vendors.rnapolis.tertiary import (
    BASE_ACCEPTORS,
    BASE_ATOMS,
    BASE_DONORS,
    BASE_EDGES,
    PHOSPHATE_ACCEPTORS,
    RIBOSE_ACCEPTORS,
    Atom,
    Residue3D,
    Structure3D,
    torsion_angle,
)

HYDROGEN_BOND_MAX_DISTANCE = 4.0
HYDROGEN_BOND_ANGLE_RANGE = (50.0, 130.0)  # 90 degrees is ideal, so allow +- 40 degrees
STACKING_MAX_DISTANCE = 6.0
STACKING_MAX_ANGLE_BETWEEN_NORMALS = 35.0
STACKING_MAX_ANGLE_BETWEEN_VECTOR_AND_NORMAL = 45.0


def angle_between_vectors(
    v1: numpy.typing.NDArray[numpy.floating],
    v2: numpy.typing.NDArray[numpy.floating],
) -> float:
    """
    Compute the angle in radians between two vectors.

    The angle is calculated using the dot product and vector norms.

    Args:
        v1: First vector.
        v2: Second vector.

    Returns:
        float: Angle between the two vectors in radians.
    """
    return math.acos(numpy.dot(v1, v2) / numpy.linalg.norm(v1) / numpy.linalg.norm(v2))


def detect_cis_trans(residue_i: Residue3D, residue_j: Residue3D) -> Optional[str]:
    """
    Classify a base pair as cis or trans based on a torsion angle.

    The function uses C1'-N1/N9 atoms of both residues to compute a torsion
    angle and classifies the pair as cis or trans.

    Args:
        residue_i: First residue in the pair.
        residue_j: Second residue in the pair.

    Returns:
        ``"c"`` for cis, ``"t"`` for trans, or ``None`` if the required atoms are missing.
    """
    c1p_i = residue_i.find_atom("C1'")
    c1p_j = residue_j.find_atom("C1'")

    if residue_i.one_letter_name in "AG":
        n9n1_i = residue_i.find_atom("N9")
    else:
        n9n1_i = residue_i.find_atom("N1")

    if residue_j.one_letter_name in "AG":
        n9n1_j = residue_j.find_atom("N9")
    else:
        n9n1_j = residue_j.find_atom("N1")

    if c1p_i is None or c1p_j is None or n9n1_i is None or n9n1_j is None:
        return None

    torsion = math.degrees(torsion_angle(c1p_i, n9n1_i, n9n1_j, c1p_j))
    return "c" if -90.0 < torsion < 90.0 else "t"


def detect_bph_br_classification(
    donor_residue: Residue3D, donor: Atom, acceptor: Atom
) -> Optional[int]:
    """
    Classify base-phosphate or base-ribose interactions.

    The classification is based on donor and acceptor atom names, residue
    type and a torsion angle around selected atoms.

    Args:
        donor_residue: Residue that contains the donor atom.
        donor: Donor atom in the interaction.
        acceptor: Acceptor atom in the interaction.

    Returns:
        Integer code describing the interaction type, or ``None`` if the interaction does not match any known class.
    """
    # source: Classification and energetics of the base-phosphate interactions
    # in RNA. Craig L. Zirbel, Judit E. Sponer, Jiri Sponer, Jesse Stombaugh
    # and Neocles B. Leontis
    if donor_residue.one_letter_name == "A":
        if donor.name == "C2":
            return 2
        if donor.name == "N6":
            n1 = donor_residue.find_atom("N1")
            c6 = donor_residue.find_atom("C6")
            if n1 is not None and c6 is not None:
                torsion = math.degrees(torsion_angle(n1, c6, donor, acceptor))
                return 6 if -90.0 < torsion < 90.0 else 7
        if donor.name == "C8":
            return 0

    if donor_residue.one_letter_name == "G":
        if donor.name == "N1":
            return 5
        if donor.name == "N2":
            n3 = donor_residue.find_atom("N3")
            c2 = donor_residue.find_atom("C2")
            if n3 is not None and c2 is not None:
                torsion = math.degrees(torsion_angle(n3, c2, donor, acceptor))
                return 1 if -90.0 < torsion < 90.0 else 3
        if donor.name == "C8":
            return 0

    if donor_residue.one_letter_name == "C":
        if donor.name == "N4":
            n3 = donor_residue.find_atom("N3")
            c4 = donor_residue.find_atom("C4")
            if n3 is not None and c4 is not None:
                torsion = math.degrees(torsion_angle(n3, c4, donor, acceptor))
                return 6 if -90.0 < torsion < 90.0 else 7
        if donor.name == "C5":
            return 9
        if donor.name == "C6":
            return 0

    if donor_residue.one_letter_name == "U":
        if donor.name == "N3":
            return 5
        if donor.name == "C5":
            return 9
        if donor.name == "C6":
            return 0

    if donor_residue.one_letter_name == "T":
        if donor.name == "N3":
            return 5
        if donor.name == "C6":
            return 0
        if donor.name == "C7":
            return 9

    return None

def find_pairs(
    structure: Structure3D, model: Optional[int] = None
) -> List[BasePair]:
    """
    Find base-base interactions.

    The function scans all residues in the structure, builds a KDTree of
    donors and acceptors and applies geometric filters to detect hydrogen
    bonds and classify them into interaction types.

    Args:
        structure: 3D structure to analyse.
        model: Optional model number to restrict the search to a single model.

    Returns:
        Lists of base pairs.
    """
    # put all donors and acceptors into a KDTree
    coordinates = []
    coordinates_atom_map: Dict[Tuple[float, float, float], Atom] = {}
    coordinates_type_map: Dict[Tuple[float, float, float], str] = {}
    coordinates_residue_map: Dict[Tuple[float, float, float], Residue3D] = {}
    for residue in structure.residues:
        if model is not None and residue.model != model:
            continue
        acceptors = (
            BASE_ACCEPTORS.get(residue.one_letter_name, [])
            + RIBOSE_ACCEPTORS
            + PHOSPHATE_ACCEPTORS
        )
        donors = BASE_DONORS.get(residue.one_letter_name, [])
        for atom_name in acceptors + donors:
            atom = residue.find_atom(atom_name)
            if atom:
                xyz = (atom.x, atom.y, atom.z)
                coordinates.append(xyz)
                coordinates_atom_map[xyz] = atom
                coordinates_type_map[xyz] = (
                    "acceptor" if atom_name in acceptors else "donor"
                )
                coordinates_residue_map[xyz] = residue

    if len(coordinates) < 2:
        return []

    kdtree = KDTree(coordinates)

    # find all hydrogen bonds
    hydrogen_bonds = []
    used_atoms: Set[Atom] = set()
    for i, j in kdtree.query_pairs(HYDROGEN_BOND_MAX_DISTANCE):
        type_i = coordinates_type_map[coordinates[i]]
        type_j = coordinates_type_map[coordinates[j]]

        # process only acceptor/donor pairs, not acceptor/acceptor or donor/donor
        if type_i == type_j:
            continue

        atom_i = coordinates_atom_map[coordinates[i]]
        atom_j = coordinates_atom_map[coordinates[j]]

        # skip spurious hydrogen bonds in the same residue
        if (
            atom_i.label is not None
            and atom_i.label is not None
            and atom_i.label == atom_j.label
        ):
            continue
        if (
            atom_i.auth is not None
            and atom_i.auth is not None
            and atom_i.auth == atom_j.auth
        ):
            continue

        residue_i = coordinates_residue_map[coordinates[i]]
        residue_j = coordinates_residue_map[coordinates[j]]

        # check for base-phosphate contacts
        ### MINIMAL VERSION: retain these checks for the update of `used_atoms`
        if (
            (atom_i.name in PHOSPHATE_ACCEPTORS or atom_j.name in PHOSPHATE_ACCEPTORS)
            and atom_i not in used_atoms
            and atom_j not in used_atoms
        ):
            if type_i == "donor":
                donor_residue = residue_i
                donor_atom, acceptor_atom = atom_i, atom_j
            else:
                donor_residue = residue_j
                donor_atom, acceptor_atom = atom_j, atom_i
            bph = detect_bph_br_classification(donor_residue, donor_atom, acceptor_atom)
            if bph is not None:
                used_atoms.add(atom_i)
                used_atoms.add(atom_j)
            continue

        # check for base-ribose contacts
        ### MINIMAL VERSION: retain these checks for the update of `used_atoms`
        if (
            (atom_i.name in RIBOSE_ACCEPTORS or atom_j.name in RIBOSE_ACCEPTORS)
            and atom_i not in used_atoms
            and atom_j not in used_atoms
        ):
            if type_i == "donor":
                donor_residue = residue_i
                donor_atom, acceptor_atom = atom_i, atom_j
            else:
                donor_residue = residue_j
                donor_atom, acceptor_atom = atom_j, atom_i
            br = detect_bph_br_classification(donor_residue, donor_atom, acceptor_atom)
            if br is not None:
                used_atoms.add(atom_i)
                used_atoms.add(atom_j)
            continue

        # check for base-base contacts
        if residue_i.base_normal_vector is None or residue_j.base_normal_vector is None:
            continue

        vector = atom_i.coordinates - atom_j.coordinates
        angle1 = math.degrees(
            angle_between_vectors(residue_i.base_normal_vector, vector)
        )
        angle2 = math.degrees(
            angle_between_vectors(residue_j.base_normal_vector, vector)
        )
        if (
            HYDROGEN_BOND_ANGLE_RANGE[0] < angle1 < HYDROGEN_BOND_ANGLE_RANGE[1]
            and HYDROGEN_BOND_ANGLE_RANGE[0] < angle2 < HYDROGEN_BOND_ANGLE_RANGE[1]
        ):
            hydrogen_bonds.append((atom_i, atom_j, residue_i, residue_j))

    # match hydrogen bonds with base edges
    labels = []
    for atom_i, atom_j, residue_i, residue_j in hydrogen_bonds:
        edges_i = BASE_EDGES.get(residue_i.one_letter_name, dict()).get(
            atom_i.name, None
        )
        edges_j = BASE_EDGES.get(residue_j.one_letter_name, dict()).get(
            atom_j.name, None
        )
        if edges_i is None or edges_j is None:
            continue

        # detect cis/trans
        cis_trans = detect_cis_trans(residue_i, residue_j)
        if cis_trans is None:
            continue

        if residue_i < residue_j:
            labels += [
                (residue_i, residue_j, cis_trans, edge_i, edge_j)
                for edge_i in edges_i for edge_j in edges_j
            ]
        else:
            labels += [
                (residue_j, residue_i, cis_trans, edge_j, edge_i)
                for edge_i in edges_i for edge_j in edges_j
            ]

    # create a list of base pairs
    base_base_pairs = []
    occupied = set()

    counter = Counter(labels)
    for interaction, hydrogen_bond_count in counter.most_common():
        if hydrogen_bond_count < 2:
            continue

        residue_i, residue_j, cis_trans, edge_i, edge_j = interaction

        if (residue_i, edge_i) in occupied:
            continue
        if (residue_j, edge_j) in occupied:
            continue

        occupied.add((residue_i, edge_i))
        occupied.add((residue_j, edge_j))

        lw = LeontisWesthof[f"{cis_trans}{edge_i}{edge_j}"]
        base_base_pairs.append((residue_i, residue_j, lw))

    base_pairs = [
        BasePair(
            Residue(residue_i.label, residue_i.auth),
            Residue(residue_j.label, residue_j.auth),
            lw,
            Saenger.from_leontis_westhof(
                residue_i.one_letter_name, residue_j.one_letter_name, lw
            ),
        )
        for residue_i, residue_j, lw in sorted(base_base_pairs)
    ]
    return base_pairs


def find_stackings(
    structure: Structure3D, model: Optional[int] = None
) -> List[Stacking]:
    """
    Find base stacking interactions in a 3D structure.

    Stacking is detected using geometric centres of bases, distance
    thresholds and angles between base normals and connecting vectors.

    Args:
        structure: 3D structure to analyse.
        model: Optional model number to restrict the search.

    Returns:
        List of detected stacking interactions.
    """
    # put all nitrogen ring centers into a KDTree
    coordinates = []
    coordinates_residue_map: Dict[Tuple[float, float, float], Residue3D] = {}
    for residue in structure.residues:
        if model is not None and residue.model != model:
            continue
        base_atoms = BASE_ATOMS.get(residue.one_letter_name, [])
        xs, ys, zs = [], [], []
        for atom_name in base_atoms:
            atom = residue.find_atom(atom_name)
            if atom is not None:
                xs.append(atom.x)
                ys.append(atom.y)
                zs.append(atom.z)
        if len(xs) > 0:
            geometric_center = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
            coordinates.append(geometric_center)
            coordinates_residue_map[geometric_center] = residue

    if len(coordinates) < 2:
        return []

    kdtree = KDTree(coordinates)

    # find all stacking interaction
    pairs = []
    for i, j in kdtree.query_pairs(STACKING_MAX_DISTANCE):
        residue_i = coordinates_residue_map[coordinates[i]]
        residue_j = coordinates_residue_map[coordinates[j]]

        # check angle between normals
        normal_i = residue_i.base_normal_vector
        normal_j = residue_j.base_normal_vector
        if normal_i is None or normal_j is None:
            continue

        angle = min(
            [
                angle_between_vectors(normal_i, normal_j),
                angle_between_vectors(-normal_i, normal_j),
            ]
        )
        if math.degrees(angle) > STACKING_MAX_ANGLE_BETWEEN_NORMALS:
            continue

        vector = numpy.array([coordinates[i][k] - coordinates[j][k] for k in (0, 1, 2)])
        angle = min(
            angle_between_vectors(vector, normal_i),
            angle_between_vectors(vector, normal_j),
        )
        if math.degrees(angle) > STACKING_MAX_ANGLE_BETWEEN_VECTOR_AND_NORMAL:
            continue

        same_direction = True if numpy.dot(normal_i, normal_j) > 0.0 else False

        if residue_i < residue_j:
            pairs.append((residue_i, residue_j, "upward" if same_direction else "inward"))
        else:
            pairs.append((residue_j, residue_i, "downward" if same_direction else "outward"))

    stackings = [
        Stacking(
            Residue(residue_i.label, residue_i.auth),
            Residue(residue_j.label, residue_j.auth),
            StackingTopology[topology]
        )
        for residue_i, residue_j, topology in sorted(pairs)
    ]
    return stackings


def extract_base_interactions(
    tertiary_structure: Structure3D, model: Optional[int] = None
) -> BaseInteractions:
    """
    Extract all base-level interactions from a 3D structure.

    This includes base pairs, base-phosphate, base-ribose and stacking
    interactions and returns them wrapped in a BaseInteractions object.

    Args:
        tertiary_structure: Input 3D structure.
        model: Optional model number to analyse.

    Returns:
        BaseInteractions: Container with all detected interactions.
    """
    base_pairs = find_pairs(tertiary_structure, model)
    stackings = find_stackings(tertiary_structure, model)
    return BaseInteractions.from_structure3d(
        tertiary_structure, base_pairs, stackings
    )
