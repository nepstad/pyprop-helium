"""
NameGenerator
=============

From a config object, construct lists of file names of a given kind.
"""

__all__ = ["GetBoundstateFilenames", "GetSingleParticleStatesFilename"]

from helium.siteconfig import SingleParticleLocations
from prefixgenerator import *
import os


def GetBoundstateFilenames(conf):
	pass


def GetSingleParticleStatesFilename(conf, model):
	"""
	Construct absolute path for file containing single 
	particle states	matching given arguments.
	"""
	#Get radial prefix
	radialPostfix = "_".join(GetRadialPostfix(conf))

	#Construct list of exisiting filenames for single-particle states
	singleStatesFiles = filter(os.path.exists, ["%s/eigenstates_sae_model_%s_%s.h5" % (loc, model, radialPostfix) for loc in SingleParticleLocations])

	#Check that we found some files
	if len(singleStatesFiles) == 0:
		raise Exception("Could not find single-particle data file!")

	return singleStatesFiles[0]

