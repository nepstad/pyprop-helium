"""
observables
===========

Extract 'observables' from wavefunction: expectation values, 
integrated/differential probabilities, etc.

"""
import copy
from numpy import real, linspace, outer, diff, zeros, double
from numpy import abs as nabs

#scipy not always available, workaround hack
try:
	from scipy.interpolate import RectBivariateSpline
except:
	RectBivariateSpline = None
	
import pyprop
from ..utils import RegisterAll, GetClassLogger
from .projectors import EigenstateProjector
from .projectors import ProductStateProjector 

@RegisterAll
class ContinuumObservables(object):
	"""
	Observables involving the total helium continuum, that is, all states
	with energies above -2.0 a.u. It is defined through the complement of
	the bound state projector, P_cont = I - P_bound

	Implements
	----------
	"""
	def __init__(self, psi, conf, ionThreshold = -2.0):
		"""Return a helium continuum observables object

		Input
		-----
		psi: a pyprop wavefunction
		conf: a pyprop config
		ionThreshold: (double) ionization energy threshold 

		"""
		self.Psi = psi
		self.Config = copy.copy(conf)
		self.IonThreshold = ionThreshold
		
		#get logger
		self.Logger = GetClassLogger(self)

		#setup bound state projector
		self.Logger.info("Setting up bound state projector...")
		self.BoundstateProjector = EigenstateProjector(conf)
		
		#keep ionization data
		self.InitialProbability = -1
		self.AbsorbedProbability = -1
		self.IonizationProbability = -1
		
		
	def Setup(self):
		"""
		"""
		#Step 1: get initial norm
		initPsi = pyprop.CreateWavefunction(self.Config)
		initPsi.Clear()
		self.Config.InitialCondition.function(initPsi, self.Config.InitialCondition)
		self.InitialProbability = initPsi.InnerProduct(initPsi).real
		
		#Step 2: calculate absorption
		self.AbsorbedProbability = self.InitialProbability - real(self.Psi.InnerProduct(self.Psi))

		#Step 3: remove projection onto bound states
		self.Logger.info("Removing bound states projection...")
		self.BoundstateProjector.RemoveProjection(self.Psi, self.IonThreshold)
		self.IonizationProbability = self.AbsorbedProbability + self.Psi.InnerProduct(self.Psi).real
		
		
	def GetIonizationProbability(self):
		"""Return ionization probability (double)
		"""
		return self.IonizationProbability 
	
	
@RegisterAll
class SingleParticleObservables(object):
	def __init__(self):
		raise NotImplementedError("Not implemented yet!")


class DoubleContinuumObservables(object):
	"""
	Observables involving the double continuum.

	Idea: buffer projections onto radial states

	Def. of double continuum should be implemented in deriving classes.

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
		#get logger
		self.Logger = GetClassLogger(self)

		#setup bound state projector
		self.Logger.info("Setting up bound state projector...")
		self.BoundstateProjector = EigenstateProjector(conf)

		#setup product state projector
		self.Logger.info("Setting up product state projector...")
		self.IsIonizedFilter = lambda E: 0.0 < E
		self.DoubleContinuumProjector = self._SetupDCProjector(conf)

		#explicitly orthogonalize on these spaces
		self.OtherProjectors = otherP
		
		#keep ionization data
		self.AbsorbedProbability = -1
		self.TotalIonizationProbability = -1
		self.DoubleIonizationProbability = -1

		#whether double ionization has been calculated already
		self.DoubleIonizationIsCalculated = False

		self.RadialProjections = None

	def _SetupDCProjector(self, conf):
		raise NotImplementedError("Please implement in derived class which specifies double continuum projectors to use.")
			
	def Setup(self):
		"""
		"""
		#Step 1: calculate absorption
		self.AbsorbedProbability = 1.0 - real(self.Psi.InnerProduct(self.Psi))

		#Step 2: remove projection onto bound states
		self.Logger.info("Removing bound states projection...")
		ionThreshold = -2.0
		self.BoundstateProjector.RemoveProjection(self.Psi, ionThreshold)
		self.TotalIonizationProbability = self.Psi.InnerProduct(self.Psi).real
		
		#Step 3: remove other projections
		for P in self.OtherProjectors:
			P.RemoveProjection(self.Psi)
			
		#Step 4: Calculate projections onto double continuum basis states
		getRadStates = \
				self.DoubleContinuumProjector.GetProjectionAllRadialStates
		self.RadialProjections = getRadStates(self.Psi)


	def GetDoubleIonizationProbability(self):
		"""Calculate double ionization probability
		"""
		#check if double ionization is already calculated, if so, 
		#just use buffered value
		if not self.DoubleIonizationIsCalculated:
			I = sum([(nabs(pop)**2).sum() for l1,l2,pop in \
				self.RadialProjections])
			self.DoubleIonizationProbability = I
			self.DoubleIonizationIsCalculated = True
			
		return self.DoubleIonizationProbability


	def GetEnergyDistribution(self, maxEnergy, numEnergyPoints):
		"""Calculate double ionization energy distribution
		"""
		#maxE = 1.0
		E = linspace(0, maxEnergy, numEnergyPoints)
		dpde = zeros((len(E), len(E)), dtype=double)

		#to store double ionization prob before interpolation
		doubleIonProb = 0

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

			#add contribution to total double ionization prob.
			doubleIonProb += sum(pop.flatten())

			#scale states with 1/dE_1 dE_2
			E1 = P.SingleStatesLeft.GetRadialEnergies(l1, self.IsIonizedFilter)
			E2 = P.SingleStatesRight.GetRadialEnergies(l2, self.IsIonizedFilter)
			pop[:-1,:-1] /= outer(diff(E1), diff(E2))
			
			#2d interpolation over all states in this shell
			interpolator = RectBivariateSpline(E1[:-1], E2[:-1], pop[:-1, :-1], kx=1, ky=1, s=0.0)

			#evaluate on given energy points, and add to total dpde
			dpde += interpolator(E, E)

		#Calculate double ionization probability to check interpolation
		absErrIonProb = abs(doubleIonProb - sum(dpde.flatten()) * diff(E)[0]**2)
		relErrIonProb = absErrIonProb/doubleIonProb
		self.Logger.debug("Integrated double ionization probability: %s" % doubleIonProb)
		if relErrIonProb > 0.01:
			warnMsg = "Integrating dP/dE1dE2 does not give correct double ionization probability"
			self.Logger.warning("%s: relerr = %s, abserr = %s." % (warnMsg, relErrIonProb, absErrIonProb))
		else :
			self.Logger.debug("Difference in double ionization probability after interpolation: %s" % absErrIonProb)

		return E, dpde


	def GetAngularDistributionCoplanar(self):
		raise NotImplementedError("Not implemented yet!")


@RegisterAll
class DoubleContinuumObservablesHePlus(DoubleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of He+ continuum states.

	"""
	def _SetupDCProjector(self, conf):
		return ProductStateProjector(conf, "he+", "he+", self.IsIonizedFilter, self.IsIonizedFilter)
	

@RegisterAll
class DoubleContinuumObservablesCoulomb175(DoubleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of Z=1.75 Coulomb waves.

	"""
	def _SetupDCProjector(self, conf):
		return ProductStateProjector(conf, "c175", "c175", self.IsIonizedFilter, self.IsIonizedFilter)
		

@RegisterAll
class DoubleContinuumPlaneWaveObservables(DoubleContinuumObservables):
	"""
	Similar to DoubleContinuumObservables, but the continuum is now defined	by
	a product of plane waves.
	"""
	def __init__(self):
		raise NotImplementedError("Not implemented yet!")

