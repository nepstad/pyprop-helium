import pyprop
from indextricks import GetCoupledIndexList
from libheliumanalysis import GetWavefunctionParticleExchange

def GetSymmetrizationIndexPairs(psi):
	"""
	Returns index pairs corresponding to symmetrization of the basis:

	L = L' L = L' and l1 = l2' and l2 = ll'

	returns a list of indices (i1, i2) 
	such that psi[i1, r1, r2] -> psi[i2, r2, r1] is particle exchange
	psi(1, 2) -> psi(2, 1)
	"""
	angularRank = 0

	coupledIndexList = GetCoupledIndexList(psi, angularRank)
	nL = len(coupledIndexList)
	
	#construct all possible pairs
	coupledPairs = [(i1, i2) for i1 in range(nL) for i2 in range(nL)]
	#filter out symmetry pairs
	symmetryFilter = lambda l, r: l.L==r.L and l.M==r.M and l.l1==r.l2 and l.l2==r.l1
	#symmetryFilter = lambda l, r: l.L==r.L and l.M==r.M and l.l1==r.l1 and l.l2==r.l2
	symmetryIndexFilter = lambda i: symmetryFilter(coupledIndexList[i[0]], coupledIndexList[i[1]])
	symmetryPairs = filter(symmetryIndexFilter, coupledPairs)

	return symmetryPairs

def GetSymmetrizedWavefunction(psi):
	"""
	Symmetrizes and anti-symmetrizes the wavefunction with respect to
	particle exchange.

	Returns a tuple of the symmetrized and anti-symmetrized wavefunction 
	(symPsi, antiSymPsi)
	"""
	pyprop.AssertSingleProc()

	sym = GetSymmetrizationIndexPairs(psi)
	exchgPsi = GetWavefunctionParticleExchange(psi, sym)

	#create symmetrized wavefunction
	symPsi = psi.Copy()
	symPsi.GetData()[:] += exchgPsi.GetData()
	symPsi.GetData()[:] *= 0.5
	
	antiSymPsi = exchgPsi
	antiSymPsi.GetData()[:] -= psi.GetData()
	antiSymPsi.GetData()[:] *= 0.5

	return symPsi, antiSymPsi

def SymmetrizeWavefunction(psi, symmetrize):
	"""
	Perform in-place symmetrization of wavefunction <psi>
	"""
	pyprop.AssertSingleProc()

	if symmetrize:
		symfactor = 1
	else:
		symfactor = -1

	sym = GetSymmetrizationIndexPairs(psi)
	exchgPsi = GetWavefunctionParticleExchange(psi, sym)

	#create symmetrized wavefunction
	exchgPsi.GetData()[:] *= symfactor
	psi.GetData()[:] += exchgPsi.GetData()
	psi.GetData()[:] *= 0.5
