import sys
sys.path.append("../../..")
import helium

sys.path.append(helium.PypropLocation)
import pyprop

import helium.core.preconditioner
import helium.core.laserfunctions
import helium.eigenvalues.eigenstates

#Update Pyprop project namespace with objects from Helium project namespace
from helium.utils import UpdatePypropProjectNamespace
UpdatePypropProjectNamespace(pyprop.ProjectNamespace)

__all__ = ["namegenerator"]
