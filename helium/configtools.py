"""
configtools
===========

"""
from numpy import max, array
from utils import RegisterAll


@RegisterAll
def Getlmax(conf):
	indexIt = conf.AngularRepresentation.index_iterator
	lmax = max(array([idx.l1 for idx in indexIt]))

	return lmax

