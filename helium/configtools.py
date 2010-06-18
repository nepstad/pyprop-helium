"""
configtools
===========

"""
import logging
from numpy import max, array
from utils import RegisterAll
import pyprop


@RegisterAll
def Getlmax(conf):
	"""Return lmax from a pyprop config instance

	Input
	-----
	conf: a pyprop config object

	Returns
	-------
	lmax: (int) maximum (single-particle) angular momentum

	"""
	indexIt = conf.AngularRepresentation.index_iterator
	lmax = max(array([idx.l1 for idx in indexIt]))

	return lmax


@RegisterAll
def GetL(conf):
	"""Return L from a pyprop config instance

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


def UpdateConfig(conf, updateParams):
	"""
	Update Pyprop config object with values given in updateParams.

	Input
	-----
	conf: pyprop config object
	updateParams: list of tuples, [('section', 'param', 'value'), ...]
	
	
	Returns
	-------
	Updated pyprop config object
	
	
	Note: Section references (i.e. through 'base') will be updated as well.

	"""
	tmpConf = pyprop.Config(conf.cfgObj)
	
	#Update config
	for section, param, val in updateParams:
		logging.info("Updating config: %s(%s): %s" % (section, param, val))
		tmpConf.SetValue(section, param, val)

	#Update config object from possible changed ConfigParser object
	newConf = pyprop.Config(tmpConf.cfgObj)
	
	return newConf
