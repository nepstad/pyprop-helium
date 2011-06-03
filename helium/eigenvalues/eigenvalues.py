"""
eigenvalues
===========

Eigenvalue calculations module. Provides high-level interface to iterative
methods for finding some eigenvalues of a large Hamiltonian

"""

import os
import sys
from numpy import where, array
import tables
import pypar
import pyprop
from pyprop.serialization import RemoveExistingDataset
from helium.namecontroller import namegenerator as NameGen
from helium.utils import RegisterAll, GetFunctionLogger

@RegisterAll
def FindEigenvaluesInverseIterationsPiram(conf):
	"""
	Calculates eigenvalues and eigenvectors for conf around a shift
	by using inverse iterations and pIRAM.

	Input
	-----
	conf: a pyprop config object

	Returns
	-------
	Piram solver object
	GMRES solver object
	Computed eigenvalues

	"""

	#Setup Problem
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#Setup shift invert solver in order to perform inverse iterations
	shiftInvertSolver = pyprop.GMRESShiftInvertSolver(prop)
	prop.Config.Arpack.inverse_iterations = True
	prop.Config.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations

	#Setup eiganvalue solver
	solver = pyprop.PiramSolver(prop)
	solver.Solve()

	#Get the converged eigenvalues
	ev = solver.Solver.GetEigenvalues().copy()
	estimates = solver.Solver.GetConvergenceEstimates().copy()
	idx = where(estimates < 0)[0]
	eigenvalues = ev[idx]

	#convert from shift inverted eigenvalues to "actual" eigenvalues
	shift = conf.Arpack.shift
	eigenvalues = 1.0 / eigenvalues + shift

	return solver, shiftInvertSolver, eigenvalues

@RegisterAll
def SaveEigenvalueSolverShiftInvert(solver, shiftInvertSolver):
	"""
	Saves the output of FindEigenvaluesNearShift, including error estimates 
	to a hdf file.
	"""
	
	logger = GetFunctionLogger()

	conf = solver.BaseProblem.Config
	L = conf.AngularRepresentation.index_iterator.L
	assert len(L) == 1
	shift = conf.Arpack.shift

	#generate filename
	filename = NameGen.GetBoundstateFilename(conf, L[0])

	#Get eigenvalue error estimates
	errorEstimatesPIRAM = solver.Solver.GetErrorEstimates()
	convergenceEstimatesEig = solver.Solver.GetConvergenceEstimates()
	errorEstimatesGMRES = shiftInvertSolver.Solver.GetErrorEstimateList()

	#Get eigenvalues
	prop = solver.BaseProblem
	E = 1.0 / array(solver.GetEigenvalues()) + shift

	#remove file if it exists
	try:
		if os.path.exists(filename):
			if pyprop.ProcId == 0:
				os.remove(filename)
	except:
		logger.error("Could not remove %s (%s)" % (filename, sys.exc_info()[1]))

	#Store eigenvalues and eigenvectors
	logger.info("Now storing eigenvectors...")
	for i in range(len(E)):
		solver.SetEigenvector(prop.psi, i)
		prop.SaveWavefunctionHDF(filename, NameGen.GetEigenvectorDatasetPath(i))

	if pyprop.ProcId == 0:
		RemoveExistingDataset(filename, "/Eig/Eigenvalues")
		RemoveExistingDataset(filename, "/Eig/ErrorEstimateListGMRES")
		RemoveExistingDataset(filename, "/Eig/ErrorEstimateListPIRAM")
		RemoveExistingDataset(filename, "/Eig/ConvergenceEstimateEig")
		h5file = tables.openFile(filename, "r+")
		try:
			#myGroup = h5file.createGroup("/", "Eig")
			myGroup = h5file.getNode("/Eig")
			h5file.createArray(myGroup, "Eigenvalues", E)
			h5file.createArray(myGroup, "ErrorEstimateListGMRES", errorEstimatesGMRES)
			h5file.createArray(myGroup, "ErrorEstimateListPIRAM", errorEstimatesPIRAM)
			h5file.createArray(myGroup, "ConvergenceEstimateEig", convergenceEstimatesEig)

			#Store config
			myGroup._v_attrs.configObject = prop.Config.cfgObj
			
			#PIRAM stats
			myGroup._v_attrs.opCount = solver.Solver.GetOperatorCount()
			myGroup._v_attrs.restartCount = solver.Solver.GetRestartCount()
			myGroup._v_attrs.orthCount = solver.Solver.GetOrthogonalizationCount()
		except:
			logger.warning("Warning: could not store eigenvalues and error estimates!")
		finally:
			h5file.close()

	pypar.barrier()

