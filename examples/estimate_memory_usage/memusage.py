# Calculate the memory usage of a Helium problem, given a Pyprop config, without
# actually setting up any potentials or wavefunctions.
# 
# NOTE: Currently it does not check whether potentials may be consolidated. For a
# reliable estimate of total memory usage, you should therefore only supply a list
# of potentials corresponding to the final list of consolidated potentials (the
# 'grid_potential_list' option under 'Propagation').
#       
#
# Example of use:
#   > import memusage
#   > import pyprop
#   > conf = pyprop.Load("config.ini")
#   > memusage.ShowProblemMemoryUsage(conf)
#

import sys
sys.path.append("../../pyprop")

import numpy 
from numpy import zeros, prod
import pyprop


def ShowProblemMemoryUsage(conf):
	"""Print the memory usage of a Pyprop problem

	"""
	#get rank info from config
	ranks = GetRanks(conf)
	angularRank = GetAngularRank(conf)
	
	#assemble problem size info
	with pyprop.EnableRedirect():
		potentialInfo = GetProblemSizeInfo(conf)
		psiShape = GetWavefunctionShape(conf)
	
	numProcs = psiShape[angularRank]
	totMem = 0
	totProcMem = 0
	for potName, potInfo in potentialInfo.iteritems():
		geomList = potInfo["geometryList"]
		rowInfo = GetPotentialRowSize(geomList)
		
		memUsage = potInfo['size'] * 16 / 1024.0**3
		totMem += memUsage
	
		allRanks = range(ranks)
		allRanks.remove(angularRank)
		procMem = max(rowInfo[angularRank]) * prod([sum(rowInfo[i]) for i in allRanks]) * 16 / 1024.0**2
		totProcMem += procMem

		pyprop.PrintOut(potName)
		pyprop.PrintOut("  Total memory usage = %.2f GB" % memUsage)
		pyprop.PrintOut("  Max proc memory usage = %.2f MB" % procMem)
	
	pyprop.PrintOut("")
	pyprop.PrintOut("Number of procs: %i" % numProcs) 
	pyprop.PrintOut("Total memory usage = %.2f GB" % totMem)
	pyprop.PrintOut("Max total proc memory usage = %.2f MB" % totProcMem)


def GetProblemSizeInfo(conf):
	"""Assemble misc. info on a Pyprop problem based on config object

	"""
	distr = pyprop.CreateDistribution(conf)
	representation = pyprop.CreateRepresentation(conf, distr)
	problemRank = representation.GetInitialShape().size
	basisList = [pyprop.CreateBasisFromRepresentation(representation.GetRepresentation(i)) for i in range(problemRank)]
	potentialInfo = {}

	for pot in conf.Propagation.grid_potential_list:
		pyprop.PrintOut("Getting info for %s" % pot) 
		configSection = conf.GetSection(pot)
		pyprop.PrintOut("  Getting geometry info...") 
		geomList = GetPotentialGeometryInfo(basisList, configSection)
		pyprop.PrintOut("  Getting shape...") 
		potShape = GetTensorPotentialShape(geomList)
		potSize = prod(potShape)
		potentialInfo[pot] = {'shape': potShape, 'size': potSize, 'geometryList': geomList}

	return potentialInfo


def GetTensorPotentialShape(geomList):
	"""Return shape of tensor potential defined by geometry list

	"""
	potentialShape = [g.GetGlobalBasisPairCount() for g in geomList]
	return potentialShape


def GetPotentialGeometryInfo(basisList, configSection):
	"""Return list of geoemetries for each potential rank

	"""
	geomList = []
	for i, basis in enumerate(basisList):
		geometryName = configSection.Get("geometry%i" % i)
		geometryInfo = basis.GetGeometryInfo(geometryName)
		geomList += [geometryInfo]

	return geomList


def GetPotentialRowSize(geomList):
	"""Get number of coupling elements for each potential rank

	"""
	rowInfo = []
	for g in geomList:
		basisPairs = g.GetBasisPairs()
		rows = zeros(g.GetGlobalBasisPairCount(), dtype=int)
		for i, j in zip(basisPairs[:,0], basisPairs[:,1]):
			rows[i] += 1

		rowInfo += [rows]

	return rowInfo


def GetAngularRank(conf):
	"""Determine angular rank of a helium config

	"""
	ranks = GetRanks(conf)
	for i in range(ranks):
		if conf.Representation.Get("representation%i" % i) == "AngularRepresentation":
			return 0
	
	raise Exception("Angular rank not found!")


def GetRanks(conf):
	"""Return number of dimension of problem defined by 'conf'

	"""
	return conf.Representation.rank


def GetWavefunctionShape(conf):
	distr = pyprop.CreateDistribution(conf)
	representation = pyprop.CreateRepresentation(conf, distr)
	return representation.GetInitialShape()


