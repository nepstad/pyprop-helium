"""
NameGenerator
=============

From a config object, construct lists of file names of a given kind.

Provides:
            GetBoundstateFilenames(conf, L)
			GetAllBoundstateFilenames(conf)
			GetSingleParticleStatesFilename(conf, model)
"""
__all__ = ["GetBoundstateFilenames", "GetAllBoundstateFilenames",
		"GetSingleParticleStatesFilename"]

from helium.siteconfig import SingleParticleLocations, BoundstateLocation
from postfixgenerator import GetRadialPostfix, GetAngularPostfix
from ..utils import RegisterAll, RegisterProjectNamespace

from pyprop import Config
from numpy import unique, int32, array
import os


@RegisterAll
def GetBoundstatesFilename(conf, L):
	"""Return filename of some angular momentum bound states.

	Generate the name of a file where bound states associated with a given
	config is located. Note that L must be specified for the filename to be
	uniquely determined.

	Input
	-----
	conf: a pyprop config object
	L: integer specifying total angular momentum

	Returns: A single file name (string)
	"""

	localConf = Config(conf.cfgObj)
	localConf.AngularRepresentation.index_iterator.L = [L]

	#Get custom postfix for filename, if specified
	customPostfix = ""
	#if hasattr(conf, "Special"):
	#	customPostfix = "_%s" % getattr(conf.Special, "custom_postfix", "")

	#Generate radial and angular postfixes from config
	radialPostfix = "_".join(GetRadialPostfix(localConf))
	angularPostfix = "_".join(GetAngularPostfix(localConf))

	#Generate bound states filename
	boundstatesFilename = "%s/boundstates_%s_%s%s.h5" % (BoundstateLocation, radialPostfix, angularPostfix, customPostfix)
	
	return boundstatesFilename


@RegisterAll
def GetAllBoundstateFilenames(conf):
	"""
	Get list of bound states files for all L's.

	Input
	-----
	conf: a pyprop config object

	Returns: a list of file names (strings)
	"""
	
	#Get unique list of Ls
	Llist = unique([L for l1, l2, L, M in conf.AngularRepresentation.index_iterator])

	#Construct file list
	fileList = [GetBoundstatesFilename(conf, L=int(L)) for L in Llist]

	#Filter out non-existing files
	filteredFileList = filter(os.path.exists, fileList)

	#Check that we found some files
	if len(filteredFileList) == 0:
		raise Exception("Could not find any bound states files!")

	return filteredFileList


@RegisterAll
def GetSingleParticleStatesFilename(conf, model):
	"""
	Construct absolute path for file containing single particle states
	matching given arguments.

	Input
	-----
	conf: a pyprop config object
	model: a string indicating the atomic model (e.g. 'h', 'he+', etc)

	Returns: file name (string)
	"""
	#Get radial prefix
	radialPostfix = "_".join(GetRadialPostfix(conf))

	#Construct list of exisiting filenames for single-particle states
	singleStatesFiles = filter(os.path.exists, ["%s/eigenstates_sae_model_%s_%s.h5" % (loc, model, radialPostfix) for loc in SingleParticleLocations])

	#Check that we found some files
	if len(singleStatesFiles) == 0:
		raise Exception("Could not find single-particle data file!")

	return singleStatesFiles[0]


def GetEigenvectorDatasetName(eigenvectorIndex):
	return "Eigenvector_%i" % eigenvectorIndex 


def GetEigenvectorDatasetPath(eigenvectorIndex):
	return "/Eig/%s" % GetEigenvectorDatasetName(eigenvectorIndex) 

