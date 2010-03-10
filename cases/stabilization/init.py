import sys
sys.path.append("../..")

import inspect
import logging
import helium

from helium.core.preconditioner import RadialTwoElectronPreconditionerIfpack
import helium.eigenvalues.eigenvalues
import helium.eigenvalues.eigenstates
import helium.core.laserfunctions
import helium.propagation.propagate
import helium.propagation.tasks
from helium.utils import UpdatePypropProjectNamespace

sys.path.append(helium.PypropLocation)
import pyprop

#import from stabilization module
import namegenerator

#Update Pyprop project namespace with objects from Helium project namespace
UpdatePypropProjectNamespace(pyprop.ProjectNamespace)

#set logging level
logging.basicConfig(level=logging.DEBUG)

