#------------------------------------------------------------------------------------
#       Eigenvalue Functions using inverse iterations with Piram or Anasazi
#------------------------------------------------------------------------------------
def FindEigenvaluesDavidson(**args):
	"""
	Finds a subset of the eigenvalues/eigenvectors using Block-Davidson
	and saves the result in outFileName
	"""
	solver, eigenvalues = FindEigenvaluesDavidsonAnasazi(**args)
	
	SaveEigenvaluesDavidson(args["outFileName"], solver, args["shift"])

	return solver


def FindEigenvaluesInverseIterations(**args):
	"""
	Finds a subset of the eigenvalues/eigenvectors using inverse iterations
	and pIRAM and saves the result in args['outFileName']
	"""
	shift = args.get("shift", -2.9)
	args.pop("shift")
	#args["shift"] = shift
	useAnasazi = args["useAnasazi"]
	outFileName = args['outFileName']

	if useAnasazi:
		solver, shiftInvertSolver, eigenvalues = FindEigenvaluesNearShiftAnasazi(shift, config=config, **args)
	else:
		solver, shiftInvertSolver, eigenvalues = FindEigenvaluesNearShift(shift, config=config, **args)
	SaveEigenvalueSolverShiftInvert(outFileName, solver, shiftInvertSolver, shift)

	#Print Statistics
	if pyprop.ProcId == 0:
		shiftInvertSolver.PrintStatistics()
	PrintOut("SHIFT = %s" % shift)
	PrintOut("GMRESERRORS = %s" % shiftInvertSolver.Solver.GetErrorEstimateList())
	PrintOut("EIGENVALUES = %s" % eigenvalues)

	return solver, shiftInvertSolver


def FindEigenvaluesNearShift(shift, **args):
	"""
	Calculates eigenvalues and eigenvectors for **args around shift
	by using inverse iterations and pIRAM
	"""

	#Setup Problem
	prop = SetupProblem(silent=True, eigenvalueShift=shift, disablePreconditioner=True, **args)

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
	eigenvalues = 1.0 / eigenvalues + shift
	PrintOut("eigenvalues = %s" % eigenvalues)

	return solver, shiftInvertSolver, eigenvalues


def GetEigenvectorDatasetName(eigenvectorIndex):
	return "Eigenvector_%i" % eigenvectorIndex 


def GetEigenvectorDatasetPath(eigenvectorIndex):
	return "/Eig/%s" % GetEigenvectorDatasetName(eigenvectorIndex) 


def SaveEigenvalueSolverShiftInvert(filename, solver, shiftInvertSolver, shift):
	"""
	Saves the output of FindEigenvaluesNearShift, including error estimates 
	to a hdf file.
	"""

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
		PrintOut("Could not remove %s (%s)" % (filename, sys.exc_info()[1]))

	#Store eigenvalues and eigenvectors
	PrintOut("Now storing eigenvectors...")
	for i in range(len(E)):
		solver.SetEigenvector(prop.psi, i)
		prop.SaveWavefunctionHDF(filename, GetEigenvectorDatasetPath(i))

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
		finally:
			h5file.close()

def SaveEigenvaluesDavidson(filename, solver, shift):
	"""
	Saves the output of FindEigenvaluesDavidson to a hdf file.
	"""

	#Get eigenvalues
	prop = solver.BaseProblem
	E = array(solver.GetEigenvalues())

	#remove file if it exists
	try:
		if os.path.exists(filename):
			if pyprop.ProcId == 0:
				os.remove(filename)
	except:
		PrintOut("Could not remove %s (%s)" % (filename, sys.exc_info()[1]))

	#Store eigenvalues and eigenvectors
	PrintOut("Now storing eigenvectors...")
	for i in range(len(E)):
		solver.SetEigenvector(prop.psi, i)
		prop.SaveWavefunctionHDF(filename, GetEigenvectorDatasetPath(i))

	if pyprop.ProcId == 0:
		RemoveExistingDataset(filename, "/Eig/Eigenvalues")
		RemoveExistingDataset(filename, "/Eig/ErrorEstimateListGMRES")
		RemoveExistingDataset(filename, "/Eig/ErrorEstimateListPIRAM")
		RemoveExistingDataset(filename, "/Eig/ConvergenceEstimateEig")
		h5file = tables.openFile(filename, "r+")
		try:
			myGroup = h5file.getNode("/Eig")
			h5file.createArray(myGroup, "Eigenvalues", E)

			#Store config
			myGroup._v_attrs.configObject = prop.Config.cfgObj
			
		finally:
			h5file.close()


#------------------------------------------------------------------------------------
#                       Eigenvalue Functions not using pyprop
#------------------------------------------------------------------------------------

def FindEigenvaluesJD(howMany, shift, tol = 1e-10, maxIter = 200, dataSetPath="/", \
	configFileName="config_eigenvalues.ini", L=0, lmax=4, outFileName = "eig_jd.h5", \
	preconType = None):
	"""
	Find some eigenvalues for a given L-subspace using Jacobi-Davidson method
	"""

	#Set up problem
	idxIt = pyprop.DefaultCoupledIndexIterator(lmax = lmax, L = L)
	prop = SetupProblem(config = configFileName, index_iterator = idxIt)

	#Set up hamilton matrix
	H_ll = SetupPotentialMatrixLL(prop,[0,1])
	H = H_ll.to_sss()

	#Set up overlap matrix
	S_ll = SetupPotentialMatrixLL(prop,[2])
	S = S_ll.to_sss()

	#Set up preconditioner
	Precon = None
	if preconType:
		Precon = preconType(H_ll)

	#Call Jacobi-Davison rountine
	numConv, E, V, numIter, numIterInner = \
		pysparse.jdsym.jdsym(H, S, Precon, howMany, shift, tol, maxIter, pysparse.itsolvers.qmrs)

	#Store eigenvalues and eigenvectors
	h5file = tables.openFile(outFileName, "w")
	try:
		myGroup = h5file.createGroup("/", "Eig")
		h5file.createArray(myGroup, "Eigenvectors", V)
		h5file.createArray(myGroup, "Eigenvalues", E)
		myGroup._v_attrs.NumberOfIterations = numIter
		myGroup._v_attrs.NumberOfInnerIterations = numIterInner
		myGroup._v_attrs.NumberOfConvergedEigs = numConv
		myGroup._v_attrs.configObject = prop.Config.cfgObj
	finally:
		h5file.close()


	return numConv, E, V, numIter, numIterInner


def FindEigenvaluesDirectDiagonalization(L=0, lmax=3, storeResult=False, checkSymmetry=False, \
	outFileName = "eig_direct.h5", config="config_eigenvalues.ini"):
	"""
	Get energies and eigenstates by direct diagonalization of L-subspace matrix
	"""
	#Coupled spherical index iterator based on given lmax and L
	index_iterator = pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=L)
	
	#Set up problem
	prop = SetupProblem(config=config, index_iterator=index_iterator)

	#Set up hamilton and overlap matrices
	HamiltonMatrix = SetupBigMatrixReal(prop, [0,1])
	OverlapMatrix = SetupBigMatrixReal(prop, [2])

	#Calculate generalized eigenvalues and eigenvectors
	print "Calculating generalized eigenvalues and eigenvectors..."
	sys.stdout.flush()
	E, V = scipy.linalg.eig(HamiltonMatrix, b=OverlapMatrix)

	#Sort eigenvalues and eigenvectors
	sortIdx = argsort(E)
	E = E[sortIdx].real
	V = V[:,sortIdx]

	#Check symmetry of eigenstates
	print "Checking symmetry of eigenstates..."
	psiShape = prop.psi.GetData().shape
	lmax = psiShape[0]
	eps = 1e-9
	symmetryList = []
	if checkSymmetry:
		for stateIdx in range(len(E)):
			
			#Check symmetry of l-component with largest norm
			lSpaceNorm = \
				[scipy.linalg.norm(V[:,stateIdx].reshape(psiShape)[i,:]) for i in range(lmax)]
			largestComponentIdx = argmax(lSpaceNorm)
			v_radial = V[:,stateIdx].reshape(psiShape)[largestComponentIdx,:,:]

			#Check if even/odd or not symmetric
			if scipy.linalg.norm(v_radial - transpose(v_radial)) < eps:
				symmetryList.append("Even")
			elif scipy.linalg.norm(v_radial + transpose(v_radial)) < eps:
				symmetryList.append("Odd")
			else:
				print "E = %s, idx = %s" % (E[stateIdx], stateIdx)
				symmetryList.append("NoSym")

	#Store result
	if storeResult:
		h5file = tables.openFile(outFileName, "w")
		try:
			myGroup = h5file.createGroup("/", "Eig")
			h5file.createArray(myGroup, "Eigenvectors", V)
			h5file.createArray(myGroup, "Eigenvalues", E)
			if symmetryList:
				h5file.createArray(myGroup, "SymmetryList", symmetryList)
			myGroup._v_attrs.configObject = prop.Config.cfgObj
		finally:
			h5file.close()

	#...or return result
	else:
		return prop, HamiltonMatrix, OverlapMatrix, E, V





#------------------------------------------------------------------------------------
#       Eigenvalue Functions using Trilinos::Anasazi
#------------------------------------------------------------------------------------

def FindEigenvaluesNearShiftAnasazi(shift, **args):
	"""
	Calculates eigenvalues and eigenvectors for **args around shift
	by using inverse iterations and Anasazi
	"""
	
	#Setup Problem
	prop = SetupProblem(silent=True, eigenvalueShift=shift, disablePreconditioner=False, **args)

	###Setup shift invert solver in order to perform inverse iterations
	prop.Config.GMRES.shift = shift
	shiftInvertSolver = pyprop.GMRESShiftInvertSolver(prop)
	prop.Config.Arpack.inverse_iterations = True
	prop.Config.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations

	solver = AnasaziSolver(prop)
	#solver.ApplyMatrix = lambda src, dst, a, b: solver.Preconditioner.Preconditioner.Multiply(src, dst)
	solver.Solve()
	
	PrintOut("GMRES Error Estimates: %s " % shiftInvertSolver.Solver.GetErrorEstimateList())
	shiftInvertSolver.PrintStatistics()

	#Get the converged eigenvalues
	eigenvalues = solver.Solver.GetEigenvalues().copy()

	#convert from shift inverted eigenvalues to "actual" eigenvalues
	eigenvalues = 1.0 / eigenvalues + shift
	PrintOut("eigenvalues = %s" % eigenvalues)

	return solver, shiftInvertSolver, eigenvalues


def FindEigenvaluesDavidsonAnasazi(shift, **args):
	"""
	Calculates eigenvalues and eigenvectors for **args around shift
	by using Block-Davidson algorithm from Anasazi
	"""
	
	#Setup Problem
	prop = SetupProblem(silent=True, eigenvalueShift=shift, disablePreconditioner=False, **args)

	#Setup preconditioner
	preconditionerName = prop.Config.Arpack.preconditioner
	Preconditioner = None
	if preconditionerName:
		preconditionerSection = prop.Config.GetSection(preconditionerName)
		Preconditioner = preconditionerSection.type(prop.psi)
		preconditionerSection.Apply(Preconditioner)
		Preconditioner.SetHamiltonianScaling(1.0)
		Preconditioner.SetOverlapScaling(-shift)
		Preconditioner.Setup(prop.Propagator)

	solver = AnasaziSolver(prop, preconditioner=Preconditioner)
	#solver.ApplyMatrix = lambda src, dst, a, b: solver.Preconditioner.Preconditioner.Multiply(src, dst)
	solver.Solve()
	
	#Get the converged eigenvalues
	eigenvalues = solver.Solver.GetEigenvalues().copy()
	PrintOut("eigenvalues = %s" % eigenvalues)

	return solver, eigenvalues




