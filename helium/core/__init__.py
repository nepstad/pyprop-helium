"""
Core
====

Core functionality for Helium calculations, such potentials.

-Potential implementations (c++)
-Time functions for time-dependent potentials (python)
-Preconditioner setup for two-electron helium problems
"""

__all__ = ["libheliumcore", "preconditioner", "laserfunctions"]

import libheliumcore
import laserfunctions
import preconditioner
