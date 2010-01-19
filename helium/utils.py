"""
Utils
=====

Misc utilities are put here for now
"""

__all__ = ["GetAngularRankIndex"]

import pyprop.core
from pyprop.core import CoupledSphericalHarmonicRepresentation as coupledSphRepr

def GetAngularRankIndex(psi):
	angIdx = ([i for i in range(psi.GetRank()) if (psi.GetRepresentation().GetRepresentation(i).__class__ == pyprop.core.CoupledSphericalHarmonicRepresentation)][0])

	#hack: convert numpy int64 to int (boost::python/numpy issue on 64bit machines)
	return int(angIdx)
