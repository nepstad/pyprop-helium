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

@RegisterAll
def GetL(conf):
	"""Return L from config file

	Looks up L from given index_interator inAngularRepresentation config
	section. If L is not unique, raise exception.

	Input
	-----
	conf: a pyprop config object

	Returns
	-------
	L: (int) total angular momentum

	Raises
	------
	If L is not uniquely specified, raise exception
	"""
	indexIterator = conf.AngularRepresentation.index_iterator
	Llist = array([idx.L for idx in indexIterator])
	L = int(Llist[0])
	if any(Llist != L):
		raise Exception("L is not unique for the given config (%s)" % Llist)
	return L



