import logging
from typing import IO, Dict, List, Optional, Tuple, Union, cast

import numpy as np
from scipy.spatial import KDTree

from volgrids._vendors.rnapolis.common import MissingResidue, ResidueAuth, ResidueLabel
from volgrids._vendors.rnapolis.tertiary import Atom, Residue3D, Structure3D

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def read_3d_structure(
    cif_or_pdb: IO[str], model: Optional[int] = None, nucleic_acid_only: bool = False
) -> Structure3D:
    """Read a PDB or mmCIF structure into a Structure3D object.

    The function auto-detects whether the input is mmCIF or PDB, parses atoms,
    groups them into residues and optionally filters to nucleic-acid residues only.

    Args:
        cif_or_pdb (IO[str]): Open file-like object with PDB or mmCIF content.
        model (Optional[int]): Model number to extract (1-based). If None, the
            first model encountered is used.
        nucleic_acid_only (bool): If True, keep only nucleic-acid residues.

    Returns:
        Structure3D: Parsed 3D structure.
    """
    atoms, modified, sequence_by_entity, is_nucleic_acid_by_entity = (
        parse_cif(cif_or_pdb) if is_cif(cif_or_pdb) else parse_pdb(cif_or_pdb)
    )
    missing_residues = (
        parse_missing_residues_cif(cif_or_pdb)
        if is_cif(cif_or_pdb)
        else parse_missing_residues_pdb(cif_or_pdb)
    )
    if not atoms:
        logger.warning("No atoms parsed from file, returning empty Structure3D.")
        return Structure3D(
            [], missing_residues=missing_residues, sequence_by_entity=sequence_by_entity
        )
    available_models = {atom.model: None for atom in atoms}
    atoms_by_model = {
        model: list(filter(lambda atom: atom.model == model, atoms))
        for model in available_models
    }
    if model is not None and model in available_models:
        atoms = atoms_by_model[model]
    else:
        atoms = atoms_by_model[list(available_models.keys())[0]]
    return group_atoms(
        atoms,
        modified,
        sequence_by_entity,
        is_nucleic_acid_by_entity,
        nucleic_acid_only,
        missing_residues,
    )


def is_cif(cif_or_pdb: IO[str]) -> bool:
    return False # [WIP] minimal rnapolis v0.16.2"


def parse_cif(
    cif: IO[str],
) -> Tuple[
    List[Atom],
    Dict[Union[ResidueLabel, ResidueAuth], str],
    Dict[str, str],
    Dict[str, bool],
]:
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")


def _read_mmcif_data(cif: IO[str]):
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")


def parse_missing_residues_cif(cif: IO[str]) -> List[MissingResidue]:
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")

def parse_missing_residues_pdb(pdb: IO[str]) -> List[MissingResidue]:
    """Parse unobserved residues reported in PDB REMARK 465 records.

    This is a best-effort parser: header and empty REMARK 465 lines are
    skipped, and only lines whose residue number parses as an integer are
    kept. Only auth-style identifiers are available in PDB files.
    """
    pdb.seek(0)

    try:
        content = pdb.read()
    except io.UnsupportedOperation:
        if hasattr(pdb, "name"):
            with open(pdb.name, "r", encoding="utf-8") as handle:
                content = handle.read()
        else:
            raise
    if isinstance(content, bytes):
        content = content.decode("utf-8")

    missing: List[MissingResidue] = []
    for line in content.splitlines():
        if not line.startswith("REMARK 465"):
            continue
        fields = line[len("REMARK 465") :].split()
        if len(fields) < 3:
            continue
        name, chain, number = fields[0], fields[1], fields[2]
        parsed_number = try_parse_int(number)
        if parsed_number is None:
            continue
        icode = fields[3] if len(fields) > 3 and fields[3] not in ("?", ".") else None
        missing.append(
            MissingResidue(
                None,
                ResidueAuth(chain, parsed_number, icode, name),
                name,
            )
        )
    return missing


def parse_pdb(
    pdb: IO[str],
) -> Tuple[
    List[Atom],
    Dict[Union[ResidueLabel, ResidueAuth], str],
    Dict[str, str],
    Dict[str, bool],
]:
    """Parse a PDB file into atoms and modification information.

    Args:
        pdb (IO[str]): Open file-like object with PDB content.

    Returns:
        tuple:
            A tuple containing:

            - **List[Atom]**: list of Atom objects.
            - **Dict[Union[ResidueLabel, ResidueAuth], str]**:
              mapping of modified residues to standard residue names.
            - **Dict[str, str]**:
              empty sequence-by-entity dictionary.
            - **Dict[str, bool]**:
              empty is-nucleic-acid-by-entity dictionary.
    """
    pdb.seek(0)
    atoms_to_process: List[Atom] = []
    modified: Dict[Union[ResidueLabel, ResidueAuth], str] = {}
    model = 1

    for line in pdb.readlines():
        if line.startswith("MODEL"):
            model = int(line[10:14].strip())
        elif line.startswith("ATOM") or line.startswith("HETATM"):
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_identifier = line[21]
            residue_number = int(line[22:26].strip())
            insertion_code = line[26] if line[26] != " " else None
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            occupancy = float(line[54:60].strip())
            auth = ResidueAuth(
                chain_identifier, residue_number, insertion_code, residue_name
            )

            atoms_to_process.append(
                Atom(None, None, auth, model, atom_name, x, y, z, occupancy)
            )
        elif line.startswith("MODRES"):
            original_name = line[12:15]
            chain_identifier = line[16]
            residue_number = int(line[18:22].strip())
            insertion_code = line[23]
            standard_residue_name = line[24:27].strip()
            auth = ResidueAuth(
                chain_identifier, residue_number, insertion_code, original_name
            )
            modified[auth] = standard_residue_name

    atoms = filter_clashing_atoms(atoms_to_process)
    return atoms, modified, {}, {}


def group_atoms(
    atoms: List[Atom],
    modified: Dict[Union[ResidueLabel, ResidueAuth], str],
    sequence_by_entity: Dict[str, str],
    is_nucleic_acid_by_entity: Dict[str, bool],
    nucleic_acid_only: bool,
    missing_residues: Optional[List[MissingResidue]] = None,
) -> Structure3D:
    """Group atoms into residues and build a Structure3D.

    This function merges atoms with the same residue identifiers, infers
    one-letter residue codes (using entity sequences and heuristics for
    nucleotides) and optionally filters to nucleic-acid residues only.

    Args:
        atoms (List[Atom]): Atoms to group.
        modified (Dict[Union[ResidueLabel, ResidueAuth], str]): Mapping for
            modified residues to their standard names.
        sequence_by_entity (Dict[str, str]): Mapping from entity ID to sequence.
        is_nucleic_acid_by_entity (Dict[str, bool]): Mapping from entity ID to
            nucleic-acid flag.
        nucleic_acid_only (bool): If True, keep only nucleic-acid residues.
        missing_residues (Optional[List[MissingResidue]]): Residues reported in
            the file header as unobserved or zero-occupancy, stored on the
            structure for downstream consumers.
    """
    if not atoms:
        return Structure3D([], missing_residues=missing_residues or [])

    key_previous = (atoms[0].label, atoms[0].auth, atoms[0].model)
    residue_atoms = [atoms[0]]
    residues: List[Residue3D] = []

    for atom in atoms[1:]:
        key = (atom.label, atom.auth, atom.model)
        if key == key_previous:
            residue_atoms.append(atom)
        else:
            label = key_previous[0]
            auth = key_previous[1]
            model = key_previous[2]
            entity_id = residue_atoms[-1].entity_id
            name = get_residue_name(auth, label, modified)
            one_letter_name = get_one_letter_name(
                entity_id, label, sequence_by_entity, name
            )

            if one_letter_name not in "ACGUTN":
                one_letter_name = detect_one_letter_name(residue_atoms)

            residues.append(
                Residue3D(
                    label, auth, model, one_letter_name, tuple(residue_atoms), name
                )
            )

            key_previous = key
            residue_atoms = [atom]

    label = key_previous[0]
    auth = key_previous[1]
    model = key_previous[2]
    entity_id = residue_atoms[-1].entity_id
    name = get_residue_name(auth, label, modified)
    one_letter_name = get_one_letter_name(entity_id, label, sequence_by_entity, name)

    if one_letter_name not in "ACGUTN":
        one_letter_name = detect_one_letter_name(residue_atoms)

    residues.append(
        Residue3D(label, auth, model, one_letter_name, tuple(residue_atoms), name)
    )

    if nucleic_acid_only:
        if is_nucleic_acid_by_entity:
            residues = [
                residue
                for residue in residues
                if (entity_id := residue.atoms[0].entity_id) is not None
                and is_nucleic_acid_by_entity.get(entity_id, False)
            ]
        else:
            residues = [residue for residue in residues if residue.is_nucleotide]

    return Structure3D(
        residues,
        missing_residues=missing_residues or [],
        sequence_by_entity=sequence_by_entity,
    )


def get_residue_name(
    auth: Optional[ResidueAuth],
    label: Optional[ResidueLabel],
    modified: Dict[Union[ResidueAuth, ResidueLabel], str],
) -> str:
    """Resolve residue name based on modification annotations and auth/label data.

    Args:
        auth (Optional[ResidueAuth]): Author-style residue identifier.
        label (Optional[ResidueLabel]): Label-style residue identifier.
        modified (Dict[Union[ResidueAuth, ResidueLabel], str]): Mapping of
            modified residues to standard names.

    Returns:
        str: Residue name (possibly standardised) or a generic placeholder.
    """
    if auth is not None and auth in modified:
        name = modified[auth].lower()
    elif label is not None and label in modified:
        name = modified[label].lower()
    elif auth is not None:
        name = auth.name
    elif label is not None:
        name = label.name
    else:
        # any nucleotide
        name = "n"
    return name


def get_one_letter_name(
    entity_id: Optional[str],
    label: Optional[ResidueLabel],
    sequence_by_entity: Dict[str, str],
    name: str,
) -> str:
    """Infer a one-letter residue code from entity sequence and residue name.

    The function first tries to read from the entity-level sequence, then falls
    back to simple heuristics based on residue name (RNA, DNA, last letter).

    Args:
        entity_id (Optional[str]): Entity ID for this residue.
        label (Optional[ResidueLabel]): Label-style residue identifier.
        sequence_by_entity (Dict[str, str]): Mapping from entity ID to sequence.
        name (str): Residue name candidate.

    Returns:
        str: One-letter residue code or a generic placeholder.
    """
    # try getting the value from _entity_poly first
    if entity_id is not None and label is not None and entity_id in sequence_by_entity:
        entity_id_str = cast(str, entity_id)
        return sequence_by_entity[entity_id_str][label.number - 1]
    # RNA
    if len(name) == 1:
        return name
    # DNA
    if len(name) == 2 and name[0].upper() == "D":
        return name[1]
    # try the last letter of the name
    if str.isalpha(name[-1]):
        return name[-1]
    # any nucleotide
    return "n"


def detect_one_letter_name(atoms: List[Atom]) -> str:
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")


def _mmcif_value(s: Optional[str]) -> Optional[str]:
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")


def try_parse_int(s: Optional[str]) -> Optional[int]:
    raise NotImplementedError("[WIP] minimal rnapolis v0.16.2")


def filter_clashing_atoms(atoms: List[Atom], clash_distance: float = 0.5) -> List[Atom]:
    """Remove duplicate and clashing atoms, keeping the higher-occupancy ones.

    The function first deduplicates atoms by (label, auth, atom name) using
    occupancy, and then uses a KDTree to detect close contacts below the given
    distance and resolves them according to occupancy.

    Args:
        atoms (List[Atom]): Input list of atoms.
        clash_distance (float): Distance threshold (in Å) for detecting clashes.

    Returns:
        Filtered list of atoms.
    """
    # First, remove duplicate atoms
    unique_atoms = {}

    for i, atom in enumerate(atoms):
        key = (atom.label, atom.auth, atom.name)
        if key not in unique_atoms or atom.occupancy > unique_atoms[key].occupancy:
            unique_atoms[key] = atom

    unique_atoms_list = list(unique_atoms.values())

    # If there are zero or one atoms, no clashes can occur
    if len(unique_atoms_list) <= 1:
        return unique_atoms_list

    # Now handle clashing atoms
    coords = np.array([(atom.x, atom.y, atom.z) for atom in unique_atoms_list])
    tree = KDTree(coords)

    pairs = tree.query_pairs(r=clash_distance)

    atoms_to_keep = set(range(len(unique_atoms_list)))

    for i, j in pairs:
        if (
            unique_atoms_list[i].occupancy is None
            or unique_atoms_list[j].occupancy is None
        ):
            continue
        if unique_atoms_list[i].occupancy > unique_atoms_list[j].occupancy:
            atoms_to_keep.discard(j)
        else:
            atoms_to_keep.discard(i)

    return [unique_atoms_list[i] for i in atoms_to_keep]
