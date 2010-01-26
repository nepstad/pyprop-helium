"""
observables
===========

Extract 'observables' from wavefunction: expectation values, 
integrated/differential probabilities, etc.

"""

from numpy import real, linspace, meshgrid, outer, diff, zeros, double
from numpy import abs as nabs
from scipy.interpolate import RectBivariateSpline
from ..utils import RegisterAll
from helium.analysis.projectors import EigenstateProjector
from helium.analysis.projectors import ProductStateProjector 

@RegisterAll
class SingleParticleObservables(object):
	pass

@RegisterAll
class DoubleContinuumObservables(object):
	"""
	Observables involving the double continuum, which is defined by a 
	product of He+ single-particle states.

	Idea: buffer projections onto radial states

	Def. of double continuum: he+ / he+

	Implements
	----------
	"""

	def __init__(self, psi, conf, otherP = []):
		"""Return a double continuum observables object

		Input
		-----
		psi: a pyprop wavefunction
		conf: a pyprop config
		otherP: list of projectors with which to orthogonalize

		"""
		self.Psi = psi
		self.BoundstateProjector = EigenstateProjector(conf)
		#isIonized = lambda E: 0.0 < E < 0.1
		self.IsIonizedFilter = lambda E: 0.0 < E < 0.1
		ionF = self.IsIonizedFilter
		self.DoubleContinuumProjector = \
				ProductStateProjector(conf, "he+", "he+", ionF, ionF)

		#explicitly orthogonalize on these spaces
		self.OtherProjectors = otherP
		
		#keep ionization data
		self.AbsorbedProbability = -1
		self.TotalIonizationProbability = -1
		self.DoubleIonizationProbability = -1

		#whether double ionization has been calculated already
		self.DoubleIonizationIsCalculated = False

		self.RadialProjections = None


	def Setup(self):
		"""
		"""
		#Step 1: calculate absorption
		self.AbsorbedProbability = 1.0 - real(self.Psi.InnerProduct(self.Psi))

		#Step 2: remove projection onto bound states
		self.BoundstateProjector.RemoveProjection(self.Psi)
		self.TotalIonizationProbability = real(self.Psi.InnerProduct(self.Psi))
		
		#Step 3: remove other projections
		for P in self.OtherProjectors:
			P.RemoveProjection(self.Psi)
			
		#Step 4: Calculate projections onto double continuum basis states
		getRadStates = \
				self.DoubleContinuumProjector.GetProjectionAllRadialStates
		self.RadialProjections = getRadStates(self.Psi)


	def GetDoubleIonizationProbability(self):
		"""
		"""

		#check if double ionization is already calculated, if so, 
		#just use buffered value
		if not self.DoubleIonizationIsCalculated:
			I = sum([(nabs(pop)**2).sum() for l1,l2,pop in \
				self.RadialProjections])
			self.DoubleIonizationProbability = I
			self.DoubleIonizationIsCalculated = True
			
		return self.DoubleIonizationProbability


	def GetEnergyDistribution(self):
		maxE = 1.0
		E = linspace(0, maxE, 300)
		dpde = zeros((len(E), len(E)), dtype=double)

		#for every l-pair (l1, l2) we have a set of (i1, i2) states
		#In order to create an approx to dp/de_1 de_2, we make a 2d 
		#interpolation for each l-shell
		#and add coherently to the dpde array
		P = self.DoubleContinuumProjector
		for l1, l2, lPop in self.RadialProjections:
			#number of states in this l-shell (matching energy filter)
			n1 = P.SingleStatesLeft.GetNumberOfStates(l1, self.IsIonizedFilter)
			n2 = P.SingleStatesRight.GetNumberOfStates(l2, self.IsIonizedFilter)

			#sum up angular momentum components
			pop = (nabs(lPop)**2).sum(axis=0).reshape(n1, n2)

			#scale states with 1/dE_1 dE_2
			E1 = P.SingleStatesLeft.GetRadialEnergies(l1, self.IsIonizedFilter)
			E2 = P.SingleStatesRight.GetRadialEnergies(l2, self.IsIonizedFilter)
			pop[:-1,:-1] /= outer(diff(E1), diff(E2))
			
			#2d interpolation over all states in this shell
			interpolator = RectBivariateSpline(E1[:-1], E2[:-1], pop[:-1, :-1], kx=1, ky=1)
			dpde += interpolator(E, E)

		#Calculate double ionization probability to check interpolation
		#lpop = [sum([p for i1, i2, p in pop]) for l1, l2, pop in populations]
		#doubleIonizationProbability = sum(lpop)
		#relErrIonProb = abs(doubleIonizationProbability - sum(dpde) * diff(E)[0]**2)/doubleIonizationProbability
		#if relErrIonProb > 0.01:
		#	warnings.warn("Integrating dP/dE1dE2 does not give correct double ionization probability: relerr = %s." % relErrIonProb, RuntimeWarning)
		#pyprop.Redirect.Disable()
		#print "Double ionization probability = ", doubleIonizationProbability
		#print "Double ionization probability (integrated) = ", sum(dpde) * diff(E)[0]**2
		#pyprop.Redirect.Enable(True)

		return E, dpde


	def GetAngularDistributionCoplanar(self):
		pass


@RegisterAll
class DoubleContinuumPlaneWaveObservables(DoubleContinuumObservables):
	"""
	Similar to DoubleContinuumObservables, but the continuum is now defined	by
	a product of plane waves.
	"""
	pass

