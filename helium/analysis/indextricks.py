"""
IndexTricks
===========
"""

from numpy import r_,  array
from helium.utils import *

def GetLocalCoupledSphericalHarmonicIndices(psi, coupledIndexFilter):
	"""
	Returns the processor local indices which corresponds to a filter on
	l1, l2, L, M
	"""
	angularRank = GetAngularRankIndex(psi)
	
	#Get info about angular representation
	repr = psi.GetRepresentation().GetRepresentation(angularRank)
	distr = psi.GetRepresentation().GetDistributedModel()
	nL = repr.GetFullShape()[0]

	#Filter global indices which are on this processor
	localStart = distr.GetLocalStartIndex(int(nL), int(angularRank))
	localEnd = localStart + nL
	isLocal = lambda i: localStart <= i < localEnd
	globalIndices = filter(isLocal, r_[:nL])

	#Filter indices where coupled index is corresponding to supplied filter
	curFilter = lambda i: coupledIndexFilter(repr.Range.GetCoupledIndex( int(i) ))
	globalFilteredIndices = filter(curFilter, globalIndices)

	#map global to local indices
	globalToLocal = lambda i: int(i - localStart)
	localFilteredIndices = map(globalToLocal, globalFilteredIndices)

	return localFilteredIndices


