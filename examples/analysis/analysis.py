"""
Example of how to use the analysis functionality in Pyprop-Helium

"""
import sys
sys.path.append("../..")
sys.path.append("../../pyprop")
import pyprop
from helium.analysis.projectors import ProductStateProjector
from helium.analysis.observables import ContinuumObservables
from helium.analysis.observables import DoubleContinuumObservables
from helium.analysis.observables import SingleContinuumObservables


#Update location of bound state files
import helium.siteconfig
import helium.namecontroller.namegenerator
helium.siteconfig.BoundstateLocation = "boundstates/"
helium.siteconfig.SingleParticleLocations += ["singleparticlestates/"]
reload(helium.namecontroller.namegenerator)


def LoadPropagatedProblem(filename):
	"""Load propagated config and wavefunction
	"""
	conf = pyprop.LoadConfigFromFile(filename)
	psi = pyprop.CreateWavefunctionFromFile(filename)

	return conf, psi


def CalculateTotalIonizationProbability(conf, psi):
	"""Calculate total ionization probability by removing projection
	of atomic bound states.

	"""

	#Set energy threshold for ionization of Helium
	ionizationTreshold = -2.0
	
	#Load propagated wavefunction and config
	#conf, psi = LoadPropagatedProblem(filename)
	
	#Setup continuum observables
	contObs = ContinuumObservables(psi, conf, ionizationTreshold)
	contObs.Setup()

	#Calculate ionization probability
	totalIonProb = contObs.GetIonizationProbability()
	absorbedProb = contObs.AbsorbedProbability

	print "\n Ionization probability = %s" % totalIonProb
	print "\n Absorbed probability = %s" % absorbedProb


def CalculateDoubleIonizationProbability(filename):
	"""Calculate double ionization probability distribution by projection onto
	He+ positive-energy eigenstates

	"""

	#Set energy threshold for ionization of Helium
	ionizationTreshold = -2.0
	
	#Load propagated wavefunction and config
	conf, psi = LoadPropagatedProblem(filename)
	
	#Setup continuum observables
	dcObs = DoubleContinuumObservablesHePlus(psi, conf, 
			ionThreshold = ionizationTreshold)
	dcObs.Setup()

	#Calculate ionization probability
	doubleIonProb = dcObs.GetDoubleIonizationProbability()

	print "\n Double ionization probability = %s" % doubleIonProb


def CalculateSingleIonizationProbability(conf, psi):
	"""Calculate single ionization probability distribution by projection onto
	He+ bound states and H continuum states

	"""

	#Set energy threshold for ionization of Helium
	ionizationTreshold = -2.0
	
	#Load propagated wavefunction and config
	#conf, psi = LoadPropagatedProblem(filename)
	
	#Setup continuum observables
	scObs = SingleContinuumObservablesHydrogenHePlus(psi, conf, 
			ionThreshold = ionizationTreshold)
	scObs.Setup()

	#Calculate ionization probability
	singleIonProb = scObs.GetSingleIonizationProbability()

	print "\n Single ionization probability = %s" % singleIonProb

	return scObs


def CalculateDoubleIonizationEnergyDistribution(conf, psi, maxEnergy=2.0,
		numPoints = 200):
	"""Calculate double ionization energy distribution by projection onto
	He+ positive-energy eigenstates

	Note: returns energy grid and energy distribution; the latter is
	numPoints x	numPoints in size

	"""

	#Set energy threshold for ionization of Helium
	ionizationTreshold = -2.0
	
	#Load propagated wavefunction and config
	#conf, psi = LoadPropagatedProblem(filename)
	
	#Setup continuum observables
	dcObs = DoubleContinuumObservablesHePlus(psi, conf, 
			ionThreshold = ionizationTreshold)
	dcObs.Setup()

	#Calculate ionization probability
	doubleIonProb = dcObs.GetDoubleIonizationProbability()
	print "\n Double ionization probability = %s" % doubleIonProb

	#Calculate energy distribution
	E, dpde = dcObs.GetEnergyDistribution(maxEnergy, numPoints)

	return E, dpde


#-----------------------------------------------------------------------------
# Implementations of specific single and double continua
#-----------------------------------------------------------------------------
class DoubleContinuumObservablesHePlus(DoubleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of He+ continuum states.

	"""
	def _SetupProjector(self, conf):
		self.IsCoulombic = True
		self.Z = 2.0
		return ProductStateProjector(conf, "HeliumPlus", "HeliumPlus", 
				self.IsIonizedFilter, self.IsIonizedFilter)


class SingleContinuumObservablesHydrogenHePlus(SingleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of He+ bound states and hydrogen continuum states.

	"""
	def _SetupProjector(self, conf):
		self.Z = 1.0
		return ProductStateProjector(conf, "Hydrogen", "HeliumPlus", 
				self.IsIonizedFilter, self.IsBoundFilter)
