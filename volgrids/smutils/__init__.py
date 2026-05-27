from ._occupancies.hbaccepts import OgHBAccepts
from ._occupancies.hbdonors import OgHBDonors
from ._occupancies.stacking import OgStacking
from ._occupancies.hydrophobic import OgHydrophobic # [TODO] hydrophilic? electrostatic?

from ._operations.app_occupancy import AppOccupancy
from ._operations.app_pwoverlap import AppPwOverlap
from ._operations.app_spheres import AppSpheres

from ._operations.residues_nucleic import ResiduesNucleic
from ._operations.chemtable_ligand import ChemTableLigand

from .app_smutils import AppSMUtils


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

OG_RADIUS_STACKING: float = 2.0
OG_RADIUS_HBA: float = 2.0
OG_RADIUS_HBD: float = 2.0
OG_RADIUS_HYDROPHOBIC: float = 2.0
# OG_RADIUS_HYDROPHILIC: float = 2.0
# OG_RADIUS_APBS: float = 2.0

DEBUG_CHEMTABLE_LIGAND: bool = False

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")
