"""
Projectors
==========

"""

__all__ = ["GetPopulationProductStates", "GetPopulationEigenstates"]

from numpy import array, int32
from indextricks import *
from libheliumanalysis import CalculatePopulationRadialProductStates


def GetPopulationProductStates(psi, singleStates1, singleStates2):
	"""
	Calculates the population of psi in a set of single electron product states

	P_i =  |< SingleState1_i(1), SingleState2_j(2) | psi(1,2) >|^2

	The projection is carried out for every combination of singlestate1 and singlestate2i
	is returned in a similar structure.

	Input
	-----
	psi: a wavefunction
	singleStates1: SingleParticleStates instance
	singleStates2: another SingleParticleStates instance

	Returns: projections onto all combinations of states, 
	"""

	#Make a copy of the wavefunction and multiply 
	#integration weights and overlap matrix
	tempPsi = psi.Copy()
	repr = psi.GetRepresentation()
	repr.MultiplyIntegrationWeights(tempPsi)
	distr = psi.GetRepresentation().GetDistributedModel()

	data = tempPsi.GetData()
	population = []

	for l1, V1 in singleStates1.IterateRadialStates():
		#print "%i/%i" % (l1, len(singleStates1))
		if V1.size == 0:
			continue

		for l2, V2 in singleStates2.IterateRadialStates():
			if V2.size == 0:
				continue

			#filter out coupled spherical harmonic indices corresponding to this l
			lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and coupledIndex.l2 == l2 
			angularIndices = array(GetLocalCoupledSphericalHarmonicIndices(psi, lfilter), dtype=int32)
			if len(angularIndices) == 0:
				continue
		
			#Get the population for every combination of v1 and v2
			projV = CalculatePopulationRadialProductStates(l1, V1, l2, V2, data, angularIndices)
			cursum = sum([p for i1, i2, p in projV])
			population.append((l1, l2, projV))

	return population


def GetPopulationEigenstates(psi, eigenstates):
	"""
	Project psi on all eigenstates and return populations
	"""
	raise NotImplementedError("Not implemented yet!")

