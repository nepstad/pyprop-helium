"""
Useful finite-difference related tools that do not require Pyprop

"""
from numpy import unique

def FindProcNumber(N, maxProcPerRank=40):
	"""
	Determine possible Nproc = q x p for given number of grid points N

	Should have N % q == 0 and N % q == 0

	"""

	possibleProcs = [p*q for p in range(1,maxProcPerRank) for q in \
		range(1,maxProcPerRank) if N%p == 0 and N%q == 0] 

	possibleProcs = unique(possibleProcs)

	return possibleProcs


