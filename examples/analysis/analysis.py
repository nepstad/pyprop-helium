"""
Example of how to use the analysis functionality in Pyprop-Helium

"""
import sys
sys.path.append("../..")
import pyprop
from helium.analysis.projectors import ProductStateProjector
from helium.analysis.observables import ContinuumObservables
from helium.analysis.observables import DoubleContinuumObservables


#Update location of bound state files
import helium.siteconfig
import helium.namecontroller
helium.siteconfig.BoundstateLocation = "boundstates/"
reload(helium.namecontroller)


def LoadPropagatedProblem(filename):
	"""Load propagated config and wavefunction
	"""
	conf = pyprop.LoadConfigFromFile(filename)
	psi = pyprop.CreateWavefunctionFromFile(filename)

	return conf, psi


def CalculateTotalIonizationProbability(filename):
	"""Calculate total ionization probability by removing projection
	of atomic bound states.

	"""

	#Set energy threshold for ionization of Helium
	ionizationTreshold = -2.0
	
	#Load propagated wavefunction and config
	conf, psi = LoadPropagatedProblem(filename)
	
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


def CalculateDoubleIonizationEnergyDistribution(filename, maxEnergy=2.0,
		numPoints = 200):
	"""Calculate double ionization energy distribution by projection onto
	He+ positive-energy eigenstates

	Note: returns energy grid and energy distribution; the latter is
	numPoints x	numPoints in size

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

	#Calculate energy distribution
	E, dpde = dcObs.GetEnergyDistribution(maxEnergy, numPoints)

	return E, dpde


#-----------------------------------------------------------------------------
# Implementations of specific double continua
#-----------------------------------------------------------------------------
class DoubleContinuumObservablesHePlus(DoubleContinuumObservables):
	"""
	Observables involving the double continuum, which is defined by a 
	product of He+ continuum states.

	"""
	def _SetupProjector(self, conf):
		self.IsCoulombic = True
		self.Z = 2.0
		return ProductStateProjector(conf, "HePlus", "HePlus", 
				self.IsIonizedFilter, self.IsIonizedFilter)
