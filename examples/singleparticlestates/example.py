import scipy.linalg
from numpy import array, r_, argsort, real, sqrt, conj, dot, zeros, double, diag
import tables

import sys
sys.path.append("../..")
from helium.configtools import UpdateConfig
from einpartikkel.utils import UpdatePypropProjectNamespace
from helium.namecontroller.postfixgenerator import GetRadialPostfix
import einpartikkel.core
from einpartikkel.core.indexiterators import FixedMLmIndexIterator
import pyprop

UpdatePypropProjectNamespace(pyprop.ProjectNamespace)


def SaveRadialEigenstates(model, configFile):
	"""Calculate and store eigenstate of hamiltionian defined by 'configFile'
	and SAE model.

	"""
	myTimers = pyprop.Timers()

	#Load config and update
	updateParams = [("SAEPotential", "base", model)]
	confIn = pyprop.Load(configFile)
	conf = UpdateConfig(confIn, updateParams)

	#Setup problem
	myTimers["Setup problem"].Start()
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	myTimers["Setup problem"].Stop()

	postfix = "_".join(GetModelPostfix(model) + GetRadialPostfix(conf))
	outputFile = "eigenstates/eigenstates_sae_%s.h5" % postfix

	#Setp eigenvalues and eigenvectors
	eigenValues, eigenVectors = SetupRadialEigenstates(prop, myTimers)
	lCount = len(eigenValues)


	#Save eigenvalues and eigenvectors to file
	myTimers["Serialization"].Start()
	if outputFile != None:
		f = tables.openFile(outputFile, "w")
		try:
			#Create main group
			eigGroup = f.createGroup(f.root, "RadialEig")

			#save config object
			eigGroup._v_attrs.configObject = prop.Config.cfgObj
			f.createArray(eigGroup, "l", r_[:lCount])

			#Save each l-shell in a separate group
			for l, (E, V) in enumerate(zip(eigenValues, eigenVectors)):
				lGroup = f.createGroup(eigGroup, "L%03i" % l)
				lGroup._v_attrs.l = l
				f.createArray(lGroup, "eigenvalues", E)
				f.createArray(lGroup, "eigenvectors", V)

		finally:
			f.close()
	myTimers["Serialization"].Stop()

	print
	print myTimers


def SetupRadialEigenstates(prop, myTimers = None):
	"""
	Finds the eigenvalues and eigenvectors of the given potentials
	of prop. From the default config file, this is the field free
	SAE He+ system.

	The eigenvalues are found by setting up a radial matrix for each l-value
	and using the generalized eigenvalue solver in scipy to find all
	eigenvalues and vectors.
	
	Eigenvalues is a list of 1-d eigenvalue arrays	
	
	"""
	assert pyprop.IsSingleProc(), "Works only on a single processor"
	
	radialRank = 1

	if not myTimers:
		myTimers = pyprop.Timers()


	eigenValues = []
	eigenVectors = []

	lCount = prop.psi.GetData().shape[0]
	potIndices = range(len(prop.Propagator.BasePropagator.PotentialList))

	#Get integration weights or overlap matrix
	myTimers["Setup overlap matrix"].Start()
	if prop.psi.GetRepresentation().IsOrthogonalBasis(radialRank):
		w = prop.psi.GetRepresentation().GetGlobalWeights(radialRank)
		S = diag(w)
	else:
		S = SetupOverlapMatrix(prop)
	myTimers["Setup overlap matrix"].Stop()

	for l in range(lCount):
		l = int(l)
		myTimers["Setup hamilton matrix"].Start()
		M = SetupRadialMatrix(prop, potIndices, l)
		myTimers["Setup hamilton matrix"].Stop()

		myTimers["Diagonalize"].Start()
		E, V = scipy.linalg.eigh(a=M, b=S)
		myTimers["Diagonalize"].Stop()

		myTimers["Postprocess eigenpairs"].Start()
		idx = argsort(real(E))
		E = real(E[idx])
		eigenValues.append(E)

		#Sort and normalize eigenvectors
		sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
		V = array([v/sNorm(v) for v in [V[:,idx[i]]
			for i in range(V.shape[1])]]).transpose()
		eigenVectors.append(V)
		myTimers["Postprocess eigenpairs"].Stop()

	return eigenValues, eigenVectors


def SetupOverlapMatrix(prop):
	"""Generate dense representation of the overlap matrix 
	"""
	overlap = prop.Propagator.BasePropagator.GeneratePotential(
			prop.Config.OverlapMatrixPotential)
	overlap.SetupStep(0.)
	matrix = SetupRadialMatrix(prop, [overlap], 0)
	return matrix


def SetupRadialMatrix(prop, whichPotentials, angularIndex):
	"""Generate a dense matrix for a fixed (l.m) value

	"""
	matrixSize = prop.psi.GetData().shape[1]
	matrix = zeros((matrixSize, matrixSize), dtype=double)

	for potNum in whichPotentials:
		if isinstance(potNum, pyprop.TensorPotential):
			potential = potNum
		else:
			potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential: %s for l = %s" % (potential.Name,
				angularIndex)

		angularBasisPairs = potential.BasisPairs[0]
		idx = [idx for idx, (i,j) in enumerate(zip(angularBasisPairs[:,0], 
			angularBasisPairs[:,1])) if i==j==angularIndex]
		if len(idx) != 1:
			raise "Invalid angular indices %s" % idx
		idx = idx[0]

		basisPairs = potential.BasisPairs[1]

		for i, (x,xp) in enumerate(basisPairs):
			indexLeft = x
			indexRight = xp
			matrix[indexLeft, indexRight] += \
				potential.PotentialData[idx, i].real

	return matrix


def GetModelPostfix(model):
	return ["model_%s" % model]
