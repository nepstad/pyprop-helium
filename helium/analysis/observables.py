"""
observables
===========

Extract 'observables' from wavefunction: expectation values, 
integrated/differential probabilities, etc.

"""
import copy
from numpy import real, linspace, outer, diff, zeros, double, array, interp, pi, linspace, unique
from numpy import int32, sqrt, exp, arctan2, cos, sin, imag, maximum, conj
from numpy import abs as nabs
import scipy.interpolate
from helium.configtools import Getlmax
from helium.analysis.indextricks import GetAngularRankIndex, GetLocalCoupledSphericalHarmonicIndices

#scipy not always available, workaround hack
try:
	from scipy.interpolate import RectBivariateSpline
except:
	RectBivariateSpline = None
	
import pyprop
from ..utils import RegisterAll, GetClassLogger
from .projectors import EigenstateProjector
from .projectors import ProductStateProjector

from above import GetCoulombPhase, AddDoubleAngularProjectionCoplanar
from tools import GetSphericalHarmonics

#------------------------------------------------------------------------------
# Continuum observables 
#------------------------------------------------------------------------------
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
		self.BoundstateProjector = EigenstateProjector(conf, ionThreshold)
		
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
		self.BoundstateProjector.RemoveProjection(self.Psi)
		self.IonizationProbability = self.AbsorbedProbability + self.Psi.InnerProduct(self.Psi).real
		
		
	def GetIonizationProbability(self):
		"""Return ionization probability (double)
		"""
		return self.IonizationProbability 
	
	
@RegisterAll
class SingleParticleObservables(object):
	def __init__(self):
		raise NotImplementedError("Not implemented yet!")

#------------------------------------------------------------------------------
# Product state observables 
#------------------------------------------------------------------------------
class ProductStateContinuumObservables(object):
	"""
	Observables involving the continuum spanned by one-electron product states.

	Def. of continuum should be implemented in deriving classes.

	Implements
	----------
	"""

	def __init__(self, psi, conf, otherP = [], ionThreshold = -2.0):
		"""Return a double continuum observables object

		Input
		-----
		psi: a pyprop wavefunction
		conf: a pyprop config
		otherP: list of projectors with which to orthogonalize

		"""
		self.Psi = psi
		self.Config = copy.copy(conf)
		self.IonizationThreshold = ionThreshold

		#get logger
		self.Logger = GetClassLogger(self)

		#setup bound state projector
		self.Logger.info("Setting up bound state projector...")
		self.BoundstateProjector = EigenstateProjector(conf, ionThreshold)

		#setup product state projector
		self.IsCoulombic = False
		self.Z = None
		self.Logger.info("Setting up product state projector...")
		self.IsIonizedFilter = lambda E: 0.0 < E
		self.IsBoundFilter = lambda E: not self.IsIonizedFilter(E)
		self.ContinuumProjector = self._SetupProjector(conf)

		#explicitly orthogonalize on these spaces
		self.OtherProjectors = otherP
		
		#keep ionization data
		self.InitialProbability = -1
		self.AbsorbedProbability = -1
		self.TotalIonizationProbability = -1
		self.DoubleIonizationProbability = -1
		self.SingleIonizationProbability = -1

		#whether ionization has been calculated already
		self.IonizationIsCalculated = False

		self.RadialProjections = None

	def _SetupProjector(self, conf):
		raise NotImplementedError("Please implement in derived class which specifies continuum projectors to use.")
			
	def Setup(self):
		"""
		"""
		#Step 1: get initial norm
		initPsi = pyprop.CreateWavefunction(self.Config)
		initPsi.Clear()
		initCond = self.Config.InitialCondition
		self.InitialProbability = 1.0
		if initCond.type == pyprop.InitialConditionType.Custom:
			initCond.function(initPsi, initCond)
			self.InitialProbability = initPsi.InnerProduct(initPsi).real
		elif initCond.type == pyprop.InitialConditionType.Function:
			localGrid = [initPsi.GetRepresentation().GetLocalGrid(i) for i in range(initPsi.GetRank())]
			initPsi.GetData()[:] = initCond.function(initCond, localGrid)
			initPsi.Normalize()
			self.InitialProbability = initPsi.InnerProduct(initPsi).real
		else:
			self.Logger.warning("Cannot handle initial condition type: %s. Cannot calculate initial norm, ionization \
			probability may be wrong." % initCond.type)
		
		#Step 2: calculate absorption
		self.AbsorbedProbability = self.InitialProbability - real(self.Psi.InnerProduct(self.Psi))

		#Step 3: remove projection onto bound states
		self.Logger.info("Removing bound states projection...")
		self.BoundstateProjector.RemoveProjection(self.Psi)
		self.TotalIonizationProbability = self.AbsorbedProbability + self.Psi.InnerProduct(self.Psi).real
		
		#Step 4: remove other projections
		for P in self.OtherProjectors:
			P.RemoveProjection(self.Psi)
			
		#Step 5: Calculate projections onto double continuum basis states
		getRadStates = self.ContinuumProjector.GetProjectionAllRadialStates
		self.RadialProjections = getRadStates(self.Psi)
		#self.RadialProjections = self.ContinuumProjector.GetPopulationProductStates(self.Psi)

#------------------------------------------------------------------------------
# Single continuum observables 
#------------------------------------------------------------------------------
class SingleContinuumObservables(ProductStateContinuumObservables):
	"""
	Observables involving the single continuum.


	Implements
	----------
	"""

	def _SetupProjector(self, conf):
		return ProductStateProjector(conf, "h", "he+", self.IsIonizedFilter, self.IsBoundFilter)
			
	def GetSingleIonizationProbability(self):
		"""Calculate single ionization probability
		"""
		#check if ionization is already calculated, if so, 
		#just use buffered value
		if not self.IonizationIsCalculated:
			#I = 2.0 * sum([sum([p for i1, i2, p in pop]) for l1, l2, pop in self.RadialProjections])
			I = 2.0 * sum([(nabs(pop)**2).sum() for l1, l2, pop in self.RadialProjections])

			self.SingleIonizationProbability = I
			self.DoubleIonizationProbability = self.TotalIonizationProbability - I
			self.IonizationIsCalculated = True

		return self.SingleIonizationProbability


	def GetEnergyDistribution(self, maxEnergy, numEnergyPoints):
		"""Calculate single ionization energy distribution
		"""

		E = linspace(0, maxEnergy, numEnergyPoints)
		dpde = zeros(len(E), dtype=double)

		#to store single ionization prob before interpolation
		singleIonProb = 0

		#for every l-pair (lIon, lBound) we have a set of (iIon, iBound) states
		#In order to create an approx to dp/de, we make an interpolation for 
		#each l-shell and add incoherently to the dpde array
		P = self.ContinuumProjector
		for lIon, lBound, lPop in self.RadialProjections:
			#number of states in this l-shell (matching energy filter)
			nIon = P.SingleStatesLeft.GetNumberOfStates(lIon, self.IsIonizedFilter)
			nBound = P.SingleStatesRight.GetNumberOfStates(lBound, self.IsBoundFilter)

			#sum up angular momentum components
			pop = 2.0 * (nabs(lPop)**2).sum(axis=0).reshape(nIon, nBound).transpose()

			#add contribution to total double ionization prob.
			singleIonProb += sum(pop.flatten())

			
			EIon = P.SingleStatesLeft.GetRadialEnergies(lIon, self.IsIonizedFilter)

			#Iterate over all bound states for this l-combination
			for iBound in range(nBound):
				#scale states with 1/diff(EIon)
				curPop = pop[iBound, :-1] / diff(EIon)
				#interpolate over ionized populations, add to total dpde
				dpde += interp(E, EIon[:-1], curPop)
				
		#Calculate single ionization probability to check interpolation
		absErrIonProb = abs(singleIonProb - sum(dpde.flatten()) * diff(E)[0])
		relErrIonProb = absErrIonProb/singleIonProb
		self.Logger.debug("Integrated single ionization probability: %s" % singleIonProb)
		if relErrIonProb > 0.01:
			warnMsg = "Integrating dP/dE does not give correct single ionization probability"
			self.Logger.warning("%s: relerr = %s, abserr = %s." % (warnMsg, relErrIonProb, absErrIonProb))
		else :
			self.Logger.debug("Difference in single ionization probability after interpolation: %s" % absErrIonProb)	
	
		return E, dpde


	def GetAngularDistributionCoplanar(self):
		raise NotImplementedError("Not implemented yet!")


#------------------------------------------------------------------------------
# Double continuum observables 
#------------------------------------------------------------------------------
class DoubleContinuumObservables(ProductStateContinuumObservables):
	"""
	Observables involving the double continuum.

	Def. of double continuum should be implemented in deriving classes.

	Implements
	----------
	"""

	def GetDoubleIonizationProbability(self):
		"""Calculate double ionization probability
		"""
		#check if double ionization is already calculated, if so, 
		#just use buffered value
		if not self.IonizationIsCalculated:
			I = sum([(nabs(pop)**2).sum() for l1,l2,pop in self.RadialProjections])
			#I = sum([sum([p for i1,i2,p in pop]) for l1,l2,pop in self.RadialProjections])
			self.DoubleIonizationProbability = I
			self.SingleIonizationProbability = self.TotalIonizationProbability - I
			self.IonizationIsCalculated = True
			
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
		P = self.ContinuumProjector
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


	def GetAngularDistributionCoplanar(self, energyGrid, thetaGrid, phi1, phi2):
		"""Calculate co-planar angular distribution.
		
		"""
		
		lmax = Getlmax(self.Config)
		angularRank = GetAngularRankIndex(self.Psi)

		#Get spherical harmonics at phi1 and phi2
		assocLegendre1 = array(GetSphericalHarmonics(lmax, thetaGrid, phi1), dtype=complex)
		assocLegendre2 = array(GetSphericalHarmonics(lmax, thetaGrid, phi2), dtype=complex)
		
		#calculate angular distr for double ionized psi, evaluated at phi1=phi2=phi
		#f = angular_distributions.GetDoubleAngularDistributionCoplanar
		#f(self.Psi, self.Z, energyGrid, self.RadialProjections, assocLegendre1, assocLegendre2, thetaGrid)
		
		#Make a copy of the wavefunction and multiply 
		#integration weights and overlap matrix
		tempPsi = self.Psi.Copy()
		repr = self.Psi.GetRepresentation()
		repr.MultiplyIntegrationWeights(tempPsi)
		angRepr = repr.GetRepresentation(angularRank)
	
		cg = pyprop.core.ClebschGordan()
	
		interpCount = len(energyGrid)
	
		thetaCount = len(thetaGrid)
		angularDistr = zeros((thetaCount, thetaCount, interpCount, interpCount), dtype=double)
	
		pop = 0
		angularDistrProj = zeros(angularDistr.shape, dtype=complex)
		P = self.ContinuumProjector
		doubleIonProb = 0
		for l1, l2, lPop in self.RadialProjections:
			#number of states in this l-shell (matching energy filter)
			n1 = P.SingleStatesLeft.GetNumberOfStates(l1, self.IsIonizedFilter)
			n2 = P.SingleStatesRight.GetNumberOfStates(l2, self.IsIonizedFilter)
	
			#sum up angular momentum components
			pop = (nabs(lPop)**2).sum(axis=0).reshape(n1, n2)
	
			#add contribution to total double ionization prob.
			doubleIonProb += sum(pop.flatten())
	
			#scale states with 1/dE_1 dE_2
			E1 = array(P.SingleStatesLeft.GetRadialEnergies(l1, self.IsIonizedFilter))
			E2 = array(P.SingleStatesRight.GetRadialEnergies(l2, self.IsIonizedFilter))
			
			#filter out coupled spherical harmonic indices. this gives us a set of L's for the given l1, l2, M
			lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and coupledIndex.l2 == l2  
			angularIndices = array(GetLocalCoupledSphericalHarmonicIndices(self.Psi, lfilter), dtype=int32)
			getIdx = lambda c: angRepr.Range.GetCoupledIndex(int(c))
			coupledIndices = map(getIdx, angularIndices)
	
			if len(angularIndices) == 0:
				continue
		
			#scale states with 1/dE_1 dE_2
			def GetDensity(curE):
				interiorSpacing = list(diff(curE)[1:])
				leftSpacing = (curE[1] - curE[0])
				rightSpacing = (curE[-1] - curE[-2])
				spacing = array([leftSpacing] + interiorSpacing + [rightSpacing])
				return 1.0 / sqrt(spacing)
			stateDensity = outer(GetDensity(E1), GetDensity(E2))
	
			#coulomb phases (-i)**(l1 + l2) * exp( sigma_l1 * sigma_l2 )
			phase1 = exp(1.0j * array([GetCoulombPhase(l1, -self.Z/curK) for curK in sqrt(2*E1)]))
			phase2 = exp(1.0j * array([GetCoulombPhase(l2, -self.Z/curK) for curK in sqrt(2*E2)]))
			phase = (-1.j)**(l1 + l2) * outer(phase1, phase2)
	
			#interpolate projection on equidistant energies and sum over L and M
			#interpProj = zeros((interpCount, interpCount), dtype=complex)
			for j in range(lPop.shape[0]):
				curRadialProj = phase * stateDensity * lPop[j,:,:]
	
				#interpolate in polar complex coordinates
				def dointerp():
					r = abs(curRadialProj)**2
					i = arctan2(imag(curRadialProj), real(curRadialProj))
					argr = cos(i)
					argi = sin(i)
					interpr = scipy.interpolate.RectBivariateSpline(E1, E2, r, kx=1, ky=1)(energyGrid, energyGrid)
					interpArgR = scipy.interpolate.RectBivariateSpline(E1, E2, argr, kx=1, ky=1)(energyGrid, energyGrid)
					interpArgI = scipy.interpolate.RectBivariateSpline(E1, E2, argi, kx=1, ky=1)(energyGrid, energyGrid)
					interpPhase = (interpArgR + 1.j*interpArgI) / sqrt(interpArgR**2 + interpArgI**2)
					curInterpProj = sqrt(maximum(interpr, 0)) * interpPhase
					return curInterpProj
				curInterpProj = dointerp()
	
				#Sum over m:
				def doSum():
					AddDoubleAngularProjectionCoplanar(angularDistrProj, assocLegendre1, assocLegendre2, curInterpProj, coupledIndices[j])
				doSum()
	
		#calculate projection for this m-shell
		angularDistr = real(angularDistrProj * conj(angularDistrProj))
	
		return angularDistr

#
# Implementations of specific double continua, such as i.e. He+ x He+ 
#
@RegisterAll
class DoubleContinuumObservablesHePlus(DoubleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of He+ continuum states.

	"""
	def _SetupProjector(self, conf):
		self.IsCoulombic = True
		self.Z = 2.0
		return ProductStateProjector(conf, "he+", "he+", self.IsIonizedFilter, self.IsIonizedFilter)

@RegisterAll
class DoubleContinuumObservablesHydrogen(DoubleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of Hydrogen continuum states.

	"""
	def _SetupProjector(self, conf):
		self.IsCoulombic = True
		self.Z = 1.0
		return ProductStateProjector(conf, "h", "h", self.IsIonizedFilter, self.IsIonizedFilter)


@RegisterAll
class DoubleContinuumObservablesCoulomb175(DoubleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of Z=1.75 Coulomb waves.

	"""
	def _SetupProjector(self, conf):
		self.IsCoulombic = True
		self.Z = 1.75
		return ProductStateProjector(conf, "c175", "c175", self.IsIonizedFilter, self.IsIonizedFilter)
	
	
@RegisterAll
class DoubleContinuumObservablesCoulombZ(DoubleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of Coulomb waves with charge Z, given in the config
	object.

	"""
	def _SetupProjector(self, conf):
		self.IsCoulombic = True
		self.Z = conf.CoulombWaves.charge
		return ProductStateProjector(conf, "c%s" % self.Z, "c%s" % self.Z, self.IsIonizedFilter, self.IsIonizedFilter)
		

@RegisterAll
class DoubleContinuumPlaneWaveObservables(DoubleContinuumObservables):
	"""
	Similar to DoubleContinuumObservables, but the continuum is now defined	by
	a product of plane waves.
	"""
	def __init__(self):
		raise NotImplementedError("Not implemented yet!")


#------------------------------------------------------------------------------
# Model state decompositions
#------------------------------------------------------------------------------
@RegisterAll
class SingleParticleStateDecomposition(object):
	"""
	"""
	def __init__(self, psi, conf, cutoff, leftModel, rightModel):
		"""Returns single particle product state decomposition

		Input
		-----
		psi: a pyprop wavefunction
		conf: a pyprop config object
		cutoff: ignore components with projection less than this
		leftModel: model to use for left state (he, he+, ...)
		rightModel: model to use for right state (he, he+, ...)

		"""
		self.Psi = psi
		self.Config = copy.copy(conf)
		self.PopulationCutoff = cutoff
		self.LeftModel = leftModel
		self.RightModel = rightModel

		self.OtherProjectors = []

		self.LeftFilter = lambda E: E < 0.0
		self.RightFilter = lambda E: E < 0.0
		
		#get logger
		self.Logger = GetClassLogger(self)

		self.Setup()

		
	def Setup(self):
		"""
		"""
		#Step 1: remove other projections
		for P in self.OtherProjectors:
			P.RemoveProjection(self.Psi)
		
		#Step 2: create projector
		self.Logger.info("Setting up product state projector...")
		self.SingleParticleProjector = \
				ProductStateProjector(self.Config, self.LeftModel, self.RightModel, \
				self.LeftFilter, self.RightFilter)

		#Step 3: Calculate projections onto product basis states
		self.Logger.info("Calculating populations...")
		#getRadStates = self.SingleParticleProjector.GetPopulationProductStates
		#self.RadialPopulations = getRadStates(self.Psi)
		getRadStates = self.SingleParticleProjector.GetProjectionAllRadialStates
		self.RadialProjections = getRadStates(self.Psi)


	def GetSingleParticleStateProjections(self):
		"""Calculate projections onto single particle product states.
		"""
		P = {}
		for l1, l2, pop in self.RadialPopulations:
			for n1,n2,p in pop:
				if p < self.PopulationCutoff:
					continue
				key = "%i,%i,%i,%i" % (n1,l1,n2,l2)
				if key in P:
					P[key] += p
				else:
					P[key] = p
			
		return P

	def GetAngularMomentumProjections(self):
		"""Calculate projections onto angular momentum product states.
		
		Sum over radial states to obtain populations in (l1,l2) product 
		states.
		
		"""
		P = {}
		for l1, l2, pop in self.RadialProjections:
			p = sum(abs(pop.flatten())**2)
			if p < self.PopulationCutoff:
				continue
			key = "%i,%i" % (l1,l2)
			if key in P:
				P[key] += p
			else:
				P[key] = p

		return P

