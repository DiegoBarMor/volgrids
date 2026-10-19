import string
from dataclasses import dataclass
from enum import Enum
from functools import cache, total_ordering
from typing import Dict, List, Optional, Tuple

class Molecule(Enum):
    """Simple classification of molecule type."""

    DNA = "DNA"
    RNA = "RNA"
    Other = "Other"


RNA_NAMES = frozenset({"A", "C", "G", "U", "I"})
DNA_NAMES = frozenset({"DA", "DC", "DG", "DT", "DU", "DI"})
SUGAR_ATOMS = frozenset({"C1'", "C2'", "C3'", "C4'", "O4'"})

class GlycosidicBond(Enum):
    """Orientation of the glycosidic bond."""

    anti = "anti"
    syn = "syn"


@total_ordering
class LeontisWesthof(Enum):
    """Leontis-Westhof base pair geometry classification."""

    cWW = "cWW"
    cWH = "cWH"
    cWS = "cWS"
    cHW = "cHW"
    cHH = "cHH"
    cHS = "cHS"
    cSW = "cSW"
    cSH = "cSH"
    cSS = "cSS"
    tWW = "tWW"
    tWH = "tWH"
    tWS = "tWS"
    tHW = "tHW"
    tHH = "tHH"
    tHS = "tHS"
    tSW = "tSW"
    tSH = "tSH"
    tSS = "tSS"

    @property
    def reverse(self):
        """Return the Leontis-Westhof class with swapped edges."""
        return LeontisWesthof[f"{self.name[0]}{self.name[2]}{self.name[1]}"]

    def __lt__(self, other):
        raise NotImplementedError("[WIP] minimal rnapolis v0.16.17")


class Saenger(Enum):
    """Saenger base pair classification."""

    I = "I"
    II = "II"
    III = "III"
    IV = "IV"
    V = "V"
    VI = "VI"
    VII = "VII"
    VIII = "VIII"
    IX = "IX"
    X = "X"
    XI = "XI"
    XII = "XII"
    XIII = "XIII"
    XIV = "XIV"
    XV = "XV"
    XVI = "XVI"
    XVII = "XVII"
    XVIII = "XVIII"
    XIX = "XIX"
    XX = "XX"
    XXI = "XXI"
    XXII = "XXII"
    XXIII = "XXIII"
    XXIV = "XXIV"
    XXV = "XXV"
    XXVI = "XXVI"
    XXVII = "XXVII"
    XXVIII = "XXVIII"

    @staticmethod
    def table() -> Dict[Tuple[str, str], str]:
        """Return mapping from (base pair, Leontis-Westhof) to Saenger class name."""
        return {
            ("AA", "tWW"): "I",
            ("AA", "tHH"): "II",
            ("GG", "tWW"): "III",
            ("GG", "tSS"): "IV",
            ("AA", "tWH"): "V",
            ("AA", "tHW"): "V",
            ("GG", "cWH"): "VI",
            ("GG", "cHW"): "VI",
            ("GG", "tWH"): "VII",
            ("GG", "tHW"): "VII",
            ("AG", "cWW"): "VIII",
            ("GA", "cWW"): "VIII",
            ("AG", "cHW"): "IX",
            ("GA", "cWH"): "IX",
            ("AG", "tWS"): "X",
            ("GA", "tSW"): "X",
            ("AG", "tHS"): "XI",
            ("GA", "tSH"): "XI",
            ("UU", "tWW"): "XII",
            ("TT", "tWW"): "XII",
            # XIII is UU/TT in tWW but donor-donor, so impossible
            # XIV and XV are both CC in tWW but donor-donor, so impossible
            ("UU", "cWW"): "XVI",
            ("TT", "cWW"): "XVI",
            ("CU", "tWW"): "XVII",
            ("UC", "tWW"): "XVII",
            ("CU", "cWW"): "XVIII",
            ("UC", "cWW"): "XVIII",
            ("CG", "cWW"): "XIX",
            ("GC", "cWW"): "XIX",
            ("AU", "cWW"): "XX",
            ("UA", "cWW"): "XX",
            ("AT", "cWW"): "XX",
            ("TA", "cWW"): "XX",
            ("AU", "tWW"): "XXI",
            ("UA", "tWW"): "XXI",
            ("AT", "tWW"): "XXI",
            ("TA", "tWW"): "XXI",
            ("CG", "tWW"): "XXII",
            ("GC", "tWW"): "XXII",
            ("AU", "cHW"): "XXIII",
            ("UA", "cWH"): "XXIII",
            ("AT", "cHW"): "XXIII",
            ("TA", "cWH"): "XXIII",
            ("AU", "tHW"): "XXIV",
            ("UA", "tWH"): "XXIV",
            ("AT", "tHW"): "XXIV",
            ("TA", "tWH"): "XXIV",
            ("AC", "tHW"): "XXV",
            ("CA", "tWH"): "XXV",
            ("AC", "tWW"): "XXVI",
            ("CA", "tWW"): "XXVI",
            ("GU", "tWW"): "XXVII",
            ("UG", "tWW"): "XXVII",
            ("GT", "tWW"): "XXVII",
            ("TG", "tWW"): "XXVII",
            ("GU", "cWW"): "XXVIII",
            ("UG", "cWW"): "XXVIII",
            ("GT", "cWW"): "XXVIII",
            ("TG", "cWW"): "XXVIII",
        }

    @classmethod
    def from_leontis_westhof(
        cls,
        residue_i_one_letter_name: str,
        residue_j_one_letter_name: str,
        lw: LeontisWesthof,
    ) -> Optional["Saenger"]:
        """Map a Leontis-Westhof class and base pair to a Saenger class.

        Args:
            residue_i_one_letter_name: One-letter code of the first base.
            residue_j_one_letter_name: One-letter code of the second base.
            lw: Leontis-Westhof classification.

        Returns:
            Matching Saenger class or ``None`` if not defined.
        """
        key = (f"{residue_i_one_letter_name}{residue_j_one_letter_name}", lw.value)
        if key in Saenger.table():
            return Saenger[Saenger.table()[key]]
        return None

    @property
    def is_canonical(self) -> bool:
        """Return True for canonical Watson-Crick or wobble pairs."""
        return self == Saenger.XIX or self == Saenger.XX or self == Saenger.XXVIII


class StackingTopology(Enum):
    """Relative orientation of stacked bases."""

    upward = "upward"
    downward = "downward"
    inward = "inward"
    outward = "outward"

class BR(Enum):
    """Base-ribose interaction classes."""

    _0 = "0BR"
    _1 = "1BR"
    _2 = "2BR"
    _3 = "3BR"
    _4 = "4BR"
    _5 = "5BR"
    _6 = "6BR"
    _7 = "7BR"
    _8 = "8BR"
    _9 = "9BR"


class BPh(Enum):
    """Base-phosphate interaction classes."""

    _0 = "0BPh"
    _1 = "1BPh"
    _2 = "2BPh"
    _3 = "3BPh"
    _4 = "4BPh"
    _5 = "5BPh"
    _6 = "6BPh"
    _7 = "7BPh"
    _8 = "8BPh"
    _9 = "9BPh"


@dataclass(frozen=True, order=True)
class ResidueLabel:
    """Label-style residue identifier (label chain/number/name)."""

    chain: Optional[str]
    number: int
    name: str


AMINO_ACID_ONE_LETTER = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass(frozen=True, order=True)
class ResidueAuth:
    """Auth-style residue identifier (auth chain/number/icode/name)."""

    chain: Optional[str]
    number: int
    icode: Optional[str]
    name: str


@dataclass(frozen=True)
@total_ordering
class Residue:
    """Unified residue identifier with both label and auth coordinates."""

    label: Optional[ResidueLabel]
    auth: Optional[ResidueAuth]

    def __lt__(self, other):
        """Compare residues by chain, number and insertion code."""
        return (self.chain or "", self.number, self.icode or " ") < (
            other.chain or "",
            other.number,
            other.icode or " ",
        )

    @property
    def chain(self) -> Optional[str]:
        """Return chain identifier from auth or label coordinates."""
        if self.auth is not None:
            return self.auth.chain
        if self.label is not None:
            return self.label.chain
        return None

    @property
    def number(self) -> Optional[int]:
        """Return residue number from auth or label coordinates."""
        if self.auth is not None:
            return self.auth.number
        if self.label is not None:
            return self.label.number
        return None

    @property
    def icode(self) -> Optional[str]:
        """Return insertion code or ``None`` if not set or blank."""
        if self.auth is not None:
            return self.auth.icode if self.auth.icode not in (" ", "?", ".") else None
        return None

    @property
    @cache
    def full_name(self) -> Optional[str]:
        """Return human-readable residue identifier (e.g. A.A/23^A)."""
        if self.auth is not None:
            if self.auth.chain is None or self.auth.chain.isspace():
                builder = f"{self.auth.name}"
            else:
                builder = f"{self.auth.chain}.{self.auth.name}"
            if len(self.auth.name) > 0 and self.auth.name[-1] in string.digits:
                builder += "/"
            builder += f"{self.auth.number}"
            if self.auth.icode:
                builder += f"^{self.auth.icode}"
            return builder
        elif self.label is not None:
            if self.label.chain is None or self.label.chain.isspace():
                builder = f"{self.label.name}"
            else:
                builder = f"{self.label.chain}.{self.label.name}"
            if len(self.label.name) > 0 and self.label.name[-1] in string.digits:
                builder += "/"
            builder += f"{self.label.number}"
            return builder
        return None


@dataclass(frozen=True, order=True)
class Interaction:
    """Base class for all interactions between two residues."""

    nt1: Residue
    nt2: Residue


@dataclass(frozen=True, order=True)
class BasePair(Interaction):
    """Base pair interaction with Leontis-Westhof and Saenger class."""

    lw: LeontisWesthof
    saenger: Optional[Saenger]


@dataclass(frozen=True, order=True)
class Stacking(Interaction):
    """Base stacking interaction."""

    topology: Optional[StackingTopology]


@dataclass(frozen=True, order=True)
class BaseRibose(Interaction):
    """Base-ribose interaction."""

    br: Optional[BR]


@dataclass(frozen=True, order=True)
class BasePhosphate(Interaction):
    """Base-phosphate interaction."""

    bph: Optional[BPh]


@dataclass(frozen=True, order=True)
class BaseInteractions:
    """Container for all base-level interactions in a structure."""

    base_pairs: List[BasePair]
    stackings: List[Stacking]
    base_ribose_interactions: List[BaseRibose]
    base_phosphate_interactions: List[BasePhosphate]

    @classmethod
    def from_structure3d(
        cls,
        structure3d: "Structure3D",
        base_pairs: List[BasePair],
        stackings: List[Stacking],
        base_ribose_interactions: List[BaseRibose],
        base_phosphate_interactions: List[BasePhosphate],
    ) -> "BaseInteractions":
        """Unify residue identifiers across interactions based on a 3D structure.

        This method ensures that all interactions share consistent auth/label
        coordinates and, when possible, fills in missing Saenger classes.

        Args:
            structure3d: 3D structure used as a reference.
            base_pairs: List of base-pair interactions.
            stackings: List of stacking interactions.
            base_ribose_interactions: List of base-ribose interactions.
            base_phosphate_interactions: List of base-phosphate interactions.

        Returns:
            BaseInteractions: Normalized interactions container.
        """
        cni2residue = {}
        cni2label = {}
        cni2auth = {}

        for residue3d in structure3d.residues:
            cni = (residue3d.chain, residue3d.number, residue3d.icode or None)
            cni2auth[cni] = residue3d.auth
            cni2label[cni] = residue3d.label
            cni2residue[cni] = residue3d

        def unify_nt(nt: Residue) -> Residue:
            if nt.auth is not None and nt.label is not None:
                return nt
            cni = (nt.chain, nt.number, nt.icode or None)
            if nt.auth is not None:
                return Residue(label=cni2label.get(cni, None), auth=nt.auth)
            if nt.label is not None:
                return Residue(label=nt.label, auth=cni2auth.get(cni, None))
            return nt

        base_pairs_new = []
        for base_pair in base_pairs:
            nt1 = unify_nt(base_pair.nt1)
            nt2 = unify_nt(base_pair.nt2)

            cni1 = (nt1.chain, nt1.number, nt1.icode or None)
            cni2 = (nt2.chain, nt2.number, nt2.icode or None)
            if cni1 not in cni2residue or cni2 not in cni2residue:
                saenger = base_pair.saenger
            else:
                saenger = base_pair.saenger or Saenger.from_leontis_westhof(
                    cni2residue[cni1].one_letter_name,
                    cni2residue[cni2].one_letter_name,
                    base_pair.lw,
                )
            if (
                nt1 != base_pair.nt1
                or nt2 != base_pair.nt2
                or saenger != base_pair.saenger
            ):
                base_pair = BasePair(nt1=nt1, nt2=nt2, lw=base_pair.lw, saenger=saenger)
            base_pairs_new.append(base_pair)

        stackings_new = []
        for stacking in stackings:
            nt1 = unify_nt(stacking.nt1)
            nt2 = unify_nt(stacking.nt2)
            if nt1 != stacking.nt1 or nt2 != stacking.nt2:
                stacking = Stacking(nt1=nt1, nt2=nt2, topology=stacking.topology)
            stackings_new.append(stacking)

        base_ribose_interactions_new = []
        for base_ribose in base_ribose_interactions:
            nt1 = unify_nt(base_ribose.nt1)
            nt2 = unify_nt(base_ribose.nt2)
            if nt1 != base_ribose.nt1 or nt2 != base_ribose.nt2:
                base_ribose = BaseRibose(nt1=nt1, nt2=nt2, br=base_ribose.br)
            base_ribose_interactions_new.append(base_ribose)

        base_phosphate_interactions_new = []
        for base_phosphate in base_phosphate_interactions:
            nt1 = unify_nt(base_phosphate.nt1)
            nt2 = unify_nt(base_phosphate.nt2)
            if nt1 != base_phosphate.nt1 or nt2 != base_phosphate.nt2:
                base_phosphate = BasePhosphate(nt1=nt1, nt2=nt2, bph=base_phosphate.bph)
            base_phosphate_interactions_new.append(base_phosphate)

        return cls(
            base_pairs=base_pairs_new,
            stackings=stackings_new,
            base_ribose_interactions=base_ribose_interactions_new,
            base_phosphate_interactions=base_phosphate_interactions_new,
        )

@dataclass
class Structure2D:
    """Secondary structure representation plus derived structural elements.

    This object collects interactions, BPSEQ, _,
    _, _ and optional inter-stem parameters.
    """

    base_pairs: List[BasePair]
    stackings: List[Stacking]
    base_ribose_interactions: List[BaseRibose]
    base_phosphate_interactions: List[BasePhosphate]
