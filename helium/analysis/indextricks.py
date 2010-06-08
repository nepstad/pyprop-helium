"""
IndexTricks
===========
"""

import sys
from numpy import r_
import pyprop
from ..utils import RegisterAll


@RegisterAll
def GetAngularRankIndex(psi):
	"""Return the angular rank index of a wavefunction
	"""
	angIdx = ([i for i in range(psi.GetRank()) if \
		(psi.GetRepresentation().GetRepresentation(i).__class__ == \
		pyprop.core.CoupledSphericalHarmonicRepresentation)][0])

	#hack: convert numpy int64 to int (boost::python/numpy issue
	#on 64bit machines)
	return int(angIdx)


@RegisterAll
def GetLocalCoupledSphericalHarmonicIndices(psi, coupledIndexFilter):
	"""
	Returns the processor local indices which corresponds to a filter on
	l1, l2, L, M

	Input
	-----
	psi: pyprop wavefunction
	coupledIndexFilter: a function defining a coupled spherical harmonics filter

	Output
	------
	A list of local filtered indices

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


@RegisterAll
def GetCoupledIndexList(psi, angularRank):
	"""
	Returns a list of the coupled indices in psi
	
	Input
	-----
	psi: pyprop wavefunction
	angularRank: (int) wavefunction angular rank
	

	Output
	------
	A list of coupled indices
	"""
	
	repr = psi.GetRepresentation().GetRepresentation(angularRank)
	distr = psi.GetRepresentation().GetDistributedModel()
	nL = repr.GetFullShape()[0]
	coupledIndexList = map(repr.Range.GetCoupledIndex, range(nL))
	
	return coupledIndexList

