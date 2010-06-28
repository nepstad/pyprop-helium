"""
Various transformations between different representations of the wavefunction.
"""

from numpy import zeros, double
import pyprop
import symmetrize
from helium import utils

@utils.RegisterAll
def ReconstructSymmetrizedRadialDensityOnGrid(datafileName, radialGrid, sym):
	"""
	Construct symmetrized 2D radial density on grid for a two-electron wavefunction
	
	Input
	-----
	datafileName: (string) name of HDF5 file where wavefunction resides (in /wavefuncion)
	radialGrid: (numpy array) 1D array of radial grid points to evaluate density on
	sym: (string) symmetry of wavefunction ('sym' or 'anti')
	
	Returns: radial density (2D numpy array of doubles)
	
	See also: ReconstructRadialDensityOnGrid
	
	"""
	#Recreate wavefunction
	psi = pyprop.CreateWavefunctionFromFile(datafileName)
	symPsi, antiSymPsi = symmetrize.GetSymmetrizedWavefunction(psi)
	
	#Determine symmetry to construct radial density for
	if sym == 'sym':
		return ReconstructRadialDensityOnGrid(symPsi, radialGrid)
	elif sym == 'anti':
		return ReconstructRadialDensityOnGrid(antiSymPsi, radialGrid)
	else:
		raise Exception("Unknown symmetry: %s! Valid options are 'sym' and 'anti'." % sym)


@utils.RegisterAll
def ReconstructRadialDensityOnGrid(psi, radialGrid):
	"""
	Construct 2D radial density on grid for a two-electron wavefunction
	
	The radial density is defined by an incoherent sum over angular
	momentum components,
	
	p(r1, r2) = \sum_(l1,l2,L,M) |f(r1,r2)|**2 r1 r2        k=(l1,l2,L,M)
	f_k(r1, r2) = \sum_(i,j) c_ijk B_i(r1)/r1 B_j(r2)/r2
	
	
	Input
	-----
	psi: pyprop wavefunction
	radialGrid: (numpy array) 1D array of radial grid points to evaluate density on
	
	Returns: radial density (2D numpy array of doubles)
	
	"""
	angularRank = 0
	radialRank = 1
	
	logger = utils.GetFunctionLogger()
	
	#Get Bspline object
	bsplineObject = psi.GetRepresentation().GetRepresentation(radialRank).GetBSplineObject()

	#Get angular index iterator
	#angularIndexIterator = conf.AngularRepresentation.index_iterator
	#angularIndexIterator = psi.GetRepresentation().GetRepresentation(angularRank)
	angularSize = psi.GetData().shape[angularRank]

	#Total number of B-splines
	numBsplines = bsplineObject.NumberOfBSplines

	#Array for final radial density
	radialDensity = zeros((radialGrid.size, radialGrid.size), dtype=double)

	#A buffer for reconstructing 
	buffer1D = zeros(radialGrid.size, dtype=complex)
	buffer2D = zeros((numBsplines, radialGrid.size), dtype=complex)
	psiSlice = zeros(numBsplines, dtype=complex)

	#Calculate radial density
	#for lIdx, st in enumerate(angularIndexIterator):
	logger.info("Constructing radial density...")
	for lIdx in range(angularSize):
		for i in range(numBsplines):
			psiSlice[:] = psi.GetData()[lIdx, i, :]
			buffer1D[:] = 0.0
			bsplineObject.ConstructFunctionFromBSplineExpansion(psiSlice, radialGrid, buffer1D)
			buffer2D[i,:] = buffer1D[:]

		for j in range(radialGrid.size):
			psiSlice[:] = buffer2D[:,j]
			buffer1D[:] = 0.0
			bsplineObject.ConstructFunctionFromBSplineExpansion(psiSlice, radialGrid, buffer1D)

			#Incoherent sum over partial waves (but coherent over b-spline coefficients)
			radialDensity[:,j] += abs(buffer1D[:])**2

	return radialDensity