#! /usr/bin/env python
import argparse
import logging
import math
import os
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy
import numpy.typing
from ordered_set import OrderedSet
from scipy.spatial import KDTree

from volgrids._vendors.rnapolis.common import (
    BR,
    BaseInteractions,
    BasePair,
    BasePhosphate,
    BaseRibose,
    BPh,
    BpSeq,
    LeontisWesthof,
    OtherInteraction,
    Residue,
    Saenger,
    Stacking,
    StackingTopology,
    Stem,
    Structure2D,
)
from volgrids._vendors.rnapolis.parser import read_3d_structure
from volgrids._vendors.rnapolis.tertiary import (
    AVERAGE_OXYGEN_PHOSPHORUS_DISTANCE_COVALENT,
    BASE_ACCEPTORS,
    BASE_ATOMS,
    BASE_DONORS,
    BASE_EDGES,
    PHOSPHATE_ACCEPTORS,
    RIBOSE_ACCEPTORS,
    Atom,
    Mapping2D3D,  # Added import
    Residue3D,
    Structure3D,
    torsion_angle,
)
from volgrids._vendors.rnapolis.util import handle_input_file

HYDROGEN_BOND_MAX_DISTANCE = 4.0
HYDROGEN_BOND_ANGLE_RANGE = (50.0, 130.0)  # 90 degrees is ideal, so allow +- 40 degrees
STACKING_MAX_DISTANCE = 6.0
STACKING_MAX_ANGLE_BETWEEN_NORMALS = 35.0
STACKING_MAX_ANGLE_BETWEEN_VECTOR_AND_NORMAL = 45.0

logging.basicConfig(level=os.getenv("LOGLEVEL", "INFO").upper())


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

    The function uses C1'–N1/N9 atoms of both residues to compute a torsion
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
    Classify base–phosphate or base–ribose interactions.

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


def merge_and_clean_bph_br(
    pairs: List[Tuple[Residue3D, Residue3D, int]],
) -> Dict[Tuple[Residue3D, Residue3D], OrderedSet[int]]:
    """
    Merge and simplify base–phosphate and base–ribose classifications.

    Multiple classifications for the same residue pair are merged into a set
    and then cleaned according to simple rules (for example 3+5 → 4).

    Args:
        pairs: List of residue pairs with their classification codes.

    Returns:
        Mapping from residue pairs to a cleaned set of classification codes.
    """
    bph_br_map: Dict[Tuple[Residue3D, Residue3D], OrderedSet[int]] = defaultdict(
        OrderedSet
    )
    for residue_i, residue_j, classification in pairs:
        bph_br_map[(residue_i, residue_j)].add(classification)
    for bphs_brs in bph_br_map.values():
        # 3BPh and 5BPh simultanously means that it is actually 4BPh
        if 3 in bphs_brs and 5 in bphs_brs:
            bphs_brs.remove(3)
            bphs_brs.remove(5)
            bphs_brs.add(4)
        # 7BPh and 9BPh simultanously means that it is actually 8BPh
        if 7 in bphs_brs and 9 in bphs_brs:
            bphs_brs.remove(7)
            bphs_brs.remove(9)
            bphs_brs.add(8)
    for key, bphs_brs in bph_br_map.items():
        if len(bphs_brs) > 1:
            bph_br_map[key] = OrderedSet([bphs_brs[0]])
    return bph_br_map


def find_pairs(
    structure: Structure3D, model: Optional[int] = None
) -> Tuple[List[BasePair], List[BasePhosphate], List[BaseRibose]]:
    """
    Find base–base, base–phosphate and base–ribose interactions.

    The function scans all residues in the structure, builds a KDTree of
    donors and acceptors and applies geometric filters to detect hydrogen
    bonds and classify them into interaction types.

    Args:
        structure: 3D structure to analyse.
        model: Optional model number to restrict the search to a single model.

    Returns:
        Lists of base pairs, base–phosphate and base–ribose interactions.
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
        return [], [], []

    kdtree = KDTree(coordinates)

    # find all hydrogen bonds
    hydrogen_bonds = []
    base_phosphate_pairs = []
    base_ribose_pairs = []
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
        logging.debug(
            f"Checking pair {residue_i.full_name} {atom_i.name} - "
            f"{residue_j.full_name} {atom_j.name}"
        )

        # check for base-phosphate contacts
        if (
            (atom_i.name in PHOSPHATE_ACCEPTORS or atom_j.name in PHOSPHATE_ACCEPTORS)
            and atom_i not in used_atoms
            and atom_j not in used_atoms
        ):
            logging.debug("Checking base-phosphate interaction")
            if type_i == "donor":
                donor_residue, acceptor_residue = residue_i, residue_j
                donor_atom, acceptor_atom = atom_i, atom_j
            else:
                donor_residue, acceptor_residue = residue_j, residue_i
                donor_atom, acceptor_atom = atom_j, atom_i
            bph = detect_bph_br_classification(donor_residue, donor_atom, acceptor_atom)
            if bph is not None:
                used_atoms.add(atom_i)
                used_atoms.add(atom_j)
                base_phosphate_pairs.append((donor_residue, acceptor_residue, bph))
            continue

        # check for base-ribose contacts
        if (
            (atom_i.name in RIBOSE_ACCEPTORS or atom_j.name in RIBOSE_ACCEPTORS)
            and atom_i not in used_atoms
            and atom_j not in used_atoms
        ):
            logging.debug("Checking base-ribose interaction")
            if type_i == "donor":
                donor_residue, acceptor_residue = residue_i, residue_j
                donor_atom, acceptor_atom = atom_i, atom_j
            else:
                donor_residue, acceptor_residue = residue_j, residue_i
                donor_atom, acceptor_atom = atom_j, atom_i
            br = detect_bph_br_classification(donor_residue, donor_atom, acceptor_atom)
            if br is not None:
                used_atoms.add(atom_i)
                used_atoms.add(atom_j)
                base_ribose_pairs.append((donor_residue, acceptor_residue, br))
            continue

        # check for base-base contacts
        if residue_i.base_normal_vector is None or residue_j.base_normal_vector is None:
            continue

        logging.debug("Checking base-base interaction")
        vector = atom_i.coordinates - atom_j.coordinates
        angle1 = math.degrees(
            angle_between_vectors(residue_i.base_normal_vector, vector)
        )
        angle2 = math.degrees(
            angle_between_vectors(residue_j.base_normal_vector, vector)
        )
        logging.debug(
            f"Angles between normals and hydrogen bond: {angle1:.2f} and {angle2:.2f}"
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

        logging.debug(
            f"Matched {residue_i.full_name} with {residue_j.full_name} as "
            f"{cis_trans} {edges_i} {edges_j}"
        )

        if residue_i < residue_j:
            for edge_i in edges_i:
                for edge_j in edges_j:
                    labels.append((residue_i, residue_j, cis_trans, edge_i, edge_j))
        else:
            for edge_i in edges_i:
                for edge_j in edges_j:
                    labels.append((residue_j, residue_i, cis_trans, edge_j, edge_i))

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

    base_pairs = []
    for residue_i, residue_j, lw in sorted(base_base_pairs):
        base_pairs.append(
            BasePair(
                Residue(residue_i.label, residue_i.auth),
                Residue(residue_j.label, residue_j.auth),
                lw,
                Saenger.from_leontis_westhof(
                    residue_i.one_letter_name, residue_j.one_letter_name, lw
                ),
            )
        )

    bph_map = merge_and_clean_bph_br(sorted(base_phosphate_pairs))
    base_phosphates = []
    for pair, bphs in bph_map.items():
        residue_i, residue_j = pair
        for bph in bphs:
            base_phosphates.append(
                BasePhosphate(
                    Residue(residue_i.label, residue_i.auth),
                    Residue(residue_j.label, residue_j.auth),
                    BPh[f"_{bph}"],
                )
            )

    br_map = merge_and_clean_bph_br(sorted(base_ribose_pairs))
    base_riboses = []
    for pair, brs in br_map.items():
        residue_i, residue_j = pair
        for br in brs:
            base_riboses.append(
                BaseRibose(
                    Residue(residue_i.label, residue_i.auth),
                    Residue(residue_j.label, residue_j.auth),
                    BR[f"_{br}"],
                )
            )

    return base_pairs, base_phosphates, base_riboses


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
            if same_direction:
                pairs.append((residue_i, residue_j, "upward"))
            else:
                pairs.append((residue_i, residue_j, "inward"))
        else:
            if same_direction:
                pairs.append((residue_j, residue_i, "downward"))
            else:
                pairs.append((residue_j, residue_i, "outward"))

    stackings = []
    for residue_i, residue_j, topology in sorted(pairs):
        nt1 = Residue(residue_i.label, residue_i.auth)
        nt2 = Residue(residue_j.label, residue_j.auth)
        stackings.append(Stacking(nt1, nt2, StackingTopology[topology]))

    return stackings


def extract_base_interactions(
    tertiary_structure: Structure3D, model: Optional[int] = None
) -> BaseInteractions:
    """
    Extract all base-level interactions from a 3D structure.

    This includes base pairs, base–phosphate, base–ribose and stacking
    interactions and returns them wrapped in a BaseInteractions object.

    Args:
        tertiary_structure: Input 3D structure.
        model: Optional model number to analyse.

    Returns:
        BaseInteractions: Container with all detected interactions.
    """
    base_pairs, base_phosphate, base_ribose = find_pairs(tertiary_structure, model)
    stackings = find_stackings(tertiary_structure, model)
    return BaseInteractions.from_structure3d(
        tertiary_structure, base_pairs, stackings, base_ribose, base_phosphate, []
    )


def generate_pymol_script(mapping: Mapping2D3D, stems: List[Stem]) -> str:
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")


def write_json(path: str, structure2d: Structure2D):
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")


def write_csv(path: str, structure2d: Structure2D):
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")

def write_bpseq(path: str, bpseq: BpSeq):
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")

def add_common_output_arguments(parser: argparse.ArgumentParser):
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")

def handle_output_arguments(
    args: argparse.Namespace,
    structure2d: Structure2D,
    mapping: Mapping2D3D,
    input_filename: str,
):
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")

def main():
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")

if __name__ == "__main__":
    main()
