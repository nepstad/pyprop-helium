execfile("small_box.py")
from pylab import *
import os



def FindBoundstatesAfterCheck(xsize,xmax,order, Ls, lmax):
	"""
	FindBoundstatesAfterCheck(xsize,xmax,Ls,lmax)

	Finds boundstates in accordance with the given parameters after
	checking if they allready exists.

	Parametres
	----------
	xsize	: xsize
	xmax	: xmax
	lmax 	: lmax
	Ls		: Ls
	order	: B-spline order
	
	"""


	locationOfBoundstates = '/home/sso059/pyprop-helium/boundstates/'


	# Check if the boundstates exists.
	modifiedLs = []
	for L in Ls:
		fileInQuestion = 'boundstates_grid_linear_xmax%i_xsize%i_order%i_angular_lmax%i_L%i_M0.h5'%(xmax,xsize,order,lmax,L)

		if not os.path.isfile(locationOfBoundstates + fileInQuestion):
			print 'The file ' + fileInQuestion + ' does not exist, and will be created.'
			modifiedLs.append(L)

 	if not len(modifiedLs) == 0:
		FindBoundstates(lmax,modifiedLs,xsize, xmax,order)
	else :
		print 'States allready exist.'
