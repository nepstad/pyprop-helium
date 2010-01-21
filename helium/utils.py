"""
Utils
=====

Misc utilities are put here for now.
"""

import sys
import pyprop.core
from pyprop.core import CoupledSphericalHarmonicRepresentation as coupledSphRepr

ProjectNamespace = []


def RegisterAll(obj):
	"""
	Function decorator for registering 'obj' in containing module's
	__all__ variable. 'obj' could be a class, or a function, as long
	as it has a __name__ and __module__ attribute.

	Based on ideas from:
	  http://groups.google.com/group/comp.lang.python/msg/11cbb03e09611b8a
	  http://code.activestate.com/recipes/576993/
	"""
	
	all = sys.modules[obj.__module__].__dict__.setdefault('__all__', [])
	if obj.__name__ not in all:  # Prevent duplicates if run from an IDE.
		all += [obj.__name__]
	return obj
RegisterAll(RegisterAll)


@RegisterAll
def RegisterProjectNamespace(obj):
	"""
	Function decorator for registering 'obj' in the ProjectNamespace
	variable for Helium, used by Pyprop to resolve config files content.
	"""
	global ProjectNamespace
	if obj.__name__ not in ProjectNamespace:  # Prevent duplicates if run from an IDE.
		ProjectNamespace += [obj]
	return obj


@RegisterAll
def UpdatePypropProjectNamespace(pypropProjNamespace):
	"""
	Update Pyprop project namespace with objects from this (Helium) package

	Input
	-----
	pypropProjNamespace: Pyprop project namespace

	Output
	------
	None
	"""
	#pypropProjNamespace += [ref for ref in ProjectNamespace if ref not in pypropProjNamespace]
	for ref in ProjectNamespace:
		if ref not in pypropProjNamespace.keys():
			pypropProjNamespace[ref.__name__] = ref


@RegisterAll
def GetAngularRankIndex(psi):
	angIdx = ([i for i in range(psi.GetRank()) if (psi.GetRepresentation().GetRepresentation(i).__class__ == pyprop.core.CoupledSphericalHarmonicRepresentation)][0])

	#hack: convert numpy int64 to int (boost::python/numpy issue on 64bit machines)
	return int(angIdx)

