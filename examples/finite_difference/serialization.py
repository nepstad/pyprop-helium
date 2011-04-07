from numpy import array, nonzero, ravel
import pyprop
from pyprop.extern.progressbar import progressbar

def LoadGroundstate(psi, conf):
	"""Load groundstate from HDF5 into psi

	The ground state radial grid should match the radial grid of psi, but
	fewer angular momentum components are permissible; only those (l1,l2,L,M)
	components of psi matching the ground state are updated.
	"""

	#Get filename
	filename = conf.filename

	#Find angular rank
	angularRank = -1
	for i in range(psi.GetRank()):
		r = psi.GetRepresentation().GetRepresentation(i)
		if isinstance(r, pyprop.core.CoupledSphericalHarmonicRepresentation):
			angularRank = i
	if angularRank < 0:
		raise Exception("No angular rank found?!")

	#Load ground state psi
	inPsi = pyprop.CreateWavefunctionFromFile(filename)

	assert (inPsi.GetData().shape[0] ==  psi.GetData().shape[0])
	assert (inPsi.GetData().shape[1] ==  psi.GetData().shape[1])

	#Get angular representation of groundstate psi
	angReprIn = inPsi.GetRepresentation().GetRepresentation(angularRank)

	#Get angular representation of "full" psi
	angReprOut = psi.GetRepresentation().GetRepresentation(angularRank)
	coupledIndexOut = array([angReprOut.Range.GetCoupledIndex(int(i)) for i in
			angReprOut.Range.GetGrid()])

	#Setup progress bar
	howMany = inPsi.GetData().shape[angularRank]
	widgets = ['Loading ground state: ', progressbar.Percentage(), ' ',
		progressbar.Bar(marker=progressbar.RotatingMarker())]
	pbar = progressbar.ProgressBar(widgets=widgets, maxval=howMany).start()

	#Iterate over angular rank and copy correct (l1,l2,L,M) portions from
	#loaded psi to "full" psi
	psi.Clear()
	for curIdx in range(howMany):
		pbar.update(curIdx)
		#Our current coupled index
		curCoupledIdxIn = angReprIn.Range.GetCoupledIndex(curIdx)

		#Matching indices in full psi
		outIdx  = nonzero(ravel(coupledIndexOut == curCoupledIdxIn))
		for idx in outIdx:
			psi.GetData()[:,:,int(idx)] = inPsi.GetData()[:,:,curIdx]

	pbar.finish()
#Put this function into Pyprop namespace
pyprop.ProjectNamespace[LoadGroundstate.__name__] = LoadGroundstate
