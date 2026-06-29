from ._occupancies.hbaccepts import OgHBAccepts
from ._occupancies.hbdonors import OgHBDonors
from ._occupancies.stacking import OgStacking
from ._occupancies.hydrophobic import OgHydrophobic # [TODO] hydrophilic? electrostatic?

from ._operations.app_occupancy import AppOccupancy
from ._operations.app_pwoverlap import AppPwOverlap
from ._operations.app_spheres import AppSpheres
from ._operations.app_boxes import AppBoxes

from ._operations.residues_nucleic import ResiduesNucleic
from ._operations.chemtable_ligand import ChemTableLigand

from .app_smutils import AppSMUtils


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

OG_STK_RADIUS: float = 2.0
OG_HBA_RADIUS: float = 2.0
OG_HBD_RADIUS: float = 2.0
# OG_APBS_RADIUS: float = 2.0
OG_HPHOB_RADIUS: float = 2.0
# OG_HPHIL_RADIUS: float = 2.0

DEBUG_CHEMTABLE_LIGAND: bool = False

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")
