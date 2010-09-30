import scipy.special
from numpy import array

def GetSphericalHarmonics(lmax, theta, phi):
	"""
	Compute spherical harmonics for every l and m up to (and including)
	a given lmax, evaluated in all (theta, phi) coordinates.
	"""
	sph = []
	for l in range(lmax+1):
		for m in range(-l, l+1):
			sph.append(scipy.special.sph_harm(m, l, phi, theta))
	return array(sph)
