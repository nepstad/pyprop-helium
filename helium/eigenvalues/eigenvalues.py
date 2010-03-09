"""
eigenvalues
===========

Eigenvalue calculations module. Provides high-level interface to iterative
methods for finding some eigenvalues of a large Hamiltonian

"""

from numpy import where
import logging
import pyprop
from pyprop import AnasaziSolver
import helium
from helium.namecontroller import namegenerator as NameGen
from helium.utils import RegisterAll

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

