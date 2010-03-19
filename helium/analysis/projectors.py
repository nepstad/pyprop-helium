"""
projectors
==========

Calculates projection of wavefunction onto various set of states.

"""
from numpy import array, int32
from indextricks import GetLocalCoupledSphericalHarmonicIndices
from helium.utils import RegisterAll
from helium.configtools import Getlmax
from helium.eigenvalues.eigenstates import Eigenstates
from singleparticle import SingleParticleStates
from libheliumanalysis import CalculatePopulationRadialProductStates
from libheliumanalysis import CalculateProjectionRadialProductStates


class Projector(object):
	"""
	Implements the action of a projection operator P: Pw -> w' and (1-P)w -> w'.
	
	In practise, this means that an input wavefunction is projected
	onto some subspace, and the result returned.
	"""
	
	def __init__(self):
		raise NotImplementedError("Please implement in derived class")

	def ProjectOnto(self, psi):
		raise NotImplementedError("Please implement in derived class")


	def ProjectOntoComplement(self, psi):
		raise NotImplementedError("Please implement in derived class")



@RegisterAll
class SymmetryProjector(Projector):
	pass


@RegisterAll
class EigenstateProjector(Projector):
	"""
	Implements a atomic eigenstate projector: P = sum_a |a><a|, where
	H|a> = E_a|a>. 
	"""

	def __init__(self, conf):
		self.Config = conf
		self.Eigenstates = Eigenstates(conf)


	def ProjectOnto(self, psi):
		pass

	def ProjectOntoComplement(self, psi):
		pass

	
	def RemoveProjection(self, psi, ionThreshold):
		"""Remove bound part of psi in-place
		"""
		prevL = -1
		
		#work buffer wavefunction
		projPsi = self.Eigenstates.GetBoundstates(0, ionThreshold)[0].Copy()
		projPsi.Clear()

		for L, E, boundPsi in self.Eigenstates.IterateBoundstates(ionThreshold):
			#Get the local indices corresponding to the local L
			LFilter = lambda idx: idx.L == L
			indexL = GetLocalCoupledSphericalHarmonicIndices(psi, LFilter)
		
			#Copy the part of psi corresponding to the current L to a 
			#single-L wavefunction to do projection.
			if not L == prevL:
				#check shape of workbuffer; create new if wrong
				if projPsi.GetData().shape != boundPsi.GetData().shape:
					projPsi = boundPsi.Copy()

				projPsi.GetData()[:] = psi.GetData()[indexL, :, :]
				prevL = L
		
			curProjList = []
			#calculate projection
			proj = projPsi.InnerProduct(boundPsi)
			curProjList.append(proj)

			#remove projection
			psi.GetData()[indexL, :, :] -= proj * boundPsi.GetData()
		

@RegisterAll
class ProductStateProjector(Projector):
	"""
	Implements a product state projector: P = sum_a |a><a|, where
	|a> = |left_a>|right_a>

	Implements
	----------
	GetPopulationProductStates(psi)
	  Calculate absolute square of projection onto all combination of
	  product states
	  
	GetProjectionRadialStates(psi, lLeft, lRight)
	  Calculate the projection onto radial angular momentum states
	"""

	def __init__(self, conf, modelLeft, modelRight, leftFilter, rightFilter):
		self.Config = conf
		self.EnergyFilterLeft = leftFilter
		self.EnergyFilterRight = rightFilter
		self.SingleStatesLeft = SingleParticleStates(modelLeft, conf)
		self.SingleStatesRight = SingleParticleStates(modelRight, conf)

		#add here: filter on L's


	def ProjectOnto(self, psi):
		raise NotImplementedError("Not implemented yet!")


	def ProjectOntoComplement(self, psi):
		raise NotImplementedError("Not implemented yet!")


	def GetPopulationProductStates(self, psi):
		"""
		Calculates the population of psi in a set of single electron product
		states:

		P_i =  |< SingleState1_i(1), SingleState2_j(2) | psi(1,2) >|^2

		The projection is carried out for every combination of singlestate1
		and singlestate2i is returned in a similar structure.

		Input
		-----
		psi: a pyprop wavefunction
		singleStates1: SingleParticleStates instance
		singleStates2: another SingleParticleStates instance

		Returns: projections onto all combinations of states

		"""
		population = []

		#Make a copy of the wavefunction and multiply 
		#integration weights and overlap matrix
		tempPsi = psi.Copy()
		repr = psi.GetRepresentation()
		repr.MultiplyIntegrationWeights(tempPsi)
		data = tempPsi.GetData()

		itLeftStates = self.SingleStatesLeft.IterateFilteredRadialStates
		itRightStates = self.SingleStatesRight.IterateFilteredRadialStates
		for l1, V1 in itLeftStates(self.EnergyFilterLeft):
			if V1.size == 0:
				continue
			
			for l2, V2 in itRightStates(self.EnergyFilterRight):
				if V2.size == 0:
					continue

				#filter out coupled spherical harmonic indices corresponding
				#to this l
				angularIndices = self.__GetFilteredAngularIndices(l1, l2, psi)

				#check that wavefunctions contained given angular momenta
				if len(angularIndices) == 0:
					continue
			
				#Get the population for every combination of v1 and v2
				calcPop = CalculatePopulationRadialProductStates
				projV = calcPop(l1, V1, l2, V2, data, angularIndices)
				#cursum = sum([p for i1, i2, p in projV])
				population.append((l1, l2, projV))

		return population


	def GetProjectionAllRadialStates(self, psi):
		"""Return (complex) radial state projections
		"""
		radialProjections = []

		#Make a copy of the wavefunction and multiply 
		#integration weights and overlap matrix
		tempPsi = psi.Copy()
		repr = psi.GetRepresentation()
		repr.MultiplyIntegrationWeights(tempPsi)
		data = tempPsi.GetData()

		#Get lmax
		lmax = Getlmax(self.Config)

		def innerLoop(radialProjections):
			#filter out coupled spherical harmonic indices corresponding
			#to this l
			angularIndices = self.__GetFilteredAngularIndices(lLeft, lRight, psi)
	
			#check that wavefunctions contained given angular momenta
			if len(angularIndices) == 0:
				return
	
			#get radial states for given angular momenta
			leftStates = self.SingleStatesLeft.GetFilteredRadialStates
			rightStates = self.SingleStatesRight.GetFilteredRadialStates
			Eleft, Vleft = leftStates(lLeft, self.EnergyFilterLeft)
			Eright, Vright = rightStates(lRight, self.EnergyFilterRight)

			if len(Eleft) == 0 or len(Eright) == 0:
				return
	
			#calculate projections. we define a function here
			#to get profiling information in python
			#print lLeft, Vleft.shape, lRight, Vright.shape, data.shape, \
			#	len(angularIndices)
			#sys.stdout.flush()
			def calculateProjectionCpp():
				return CalculateProjectionRadialProductStates(lLeft, Vleft, \
						lRight, Vright, data, angularIndices)
			projV = calculateProjectionCpp()
			
			radialProjections += [(lLeft, lRight, projV)]


		#iterate over all lLeft, lRight combinations
		for lLeft in range(lmax+1):
			for lRight in range(lmax+1):
				
				#calculate radial projections
				innerLoop(radialProjections)

		return radialProjections


	def __GetFilteredAngularIndices(self, l1, l2, psi):
		#filter out coupled spherical harmonic indices corresponding
		#to this l
		lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and \
			coupledIndex.l2 == l2 
		angularIndices = \
			GetLocalCoupledSphericalHarmonicIndices(psi, lfilter)
		angularIndices = array(angularIndices, dtype = int32)

		return angularIndices

