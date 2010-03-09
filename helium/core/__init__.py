"""
Core
====

Core functionality for Helium calculations, such potentials.

-Potential implementations (c++)
-Time functions for time-dependent potentials (python)
-Preconditioner setup for two-electron helium problems
"""

#put boost::python-wrapped potentials into project namespace
from helium.utils import ProjectNamespace, RegisterProjectNamespace
import libheliumcore
from libheliumcore import *
for key in libheliumcore.__dict__.iterkeys():
	if not key.startswith("__"):
		RegisterProjectNamespace(eval(key))

__all__ = ["libheliumcore", "preconditioner", "laserfunctions"]
