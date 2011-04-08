"""
PostfixGenerator
===============

From a config object, construct set of parameter prefixes to uniquely
determine some aspect of the config file in the context of Helium.

Provides:

	GetProjectPostfix(conf)
	GetRadialPostfix(conf)
	GetAngularPostifx(conf)
"""

from ..utils import RegisterAll
from numpy import unique, diff
import pyprop

@RegisterAll
def GetProjectPostfix(conf):
	"""
	Get project prefix from config file
	"""
	return conf.Description.project


@RegisterAll
def GetRadialPostfix(conf):
	"""
	Returns a "unique" list of strings string identifying the radial grid
	implied by the specified args
	"""
	cfg = conf.RadialRepresentation
	repType = cfg.type()
	if isinstance(repType, pyprop.core.BSplineRepresentation):
		return GetRadialPostfixBsplines(conf)
	elif isinstance(repType, pyprop.core.CustomGridRepresentation):
		return GetRadialPostfixFiniteDifference(conf)
	else:
		raise Exception("Unknown representation: %s" % repType)


def GetRadialPostfixBsplines(conf):
	"""Bspline radial grid postfix

	Returns a "unique" list of strings string identifying the B-spline
	radial grid	implied by the specified args
	"""
	cfg = conf.RadialRepresentation
	gridType = cfg.bpstype
	postfix = ["bsplines", "grid", gridType, "xmax%i" % cfg.xmax, 
			"xsize%i" % cfg.xsize, "order%i" % cfg.order]
	if gridType == "linear":
		pass
	elif gridType == "exponentiallinear":
		postfix.append("xpartition%i" % cfg.xpartition)
		postfix.append("gamma%.1f" % cfg.gamma)
	elif gridType == "exponential":
		postfix.append("gamma%.1f" % cfg.gamma)

	return postfix


def GetRadialPostfixFiniteDifference(conf):
	"""Finite difference radial grid postfix

	Returns a "unique" list of strings string identifying the finite difference
	radial grid	implied by the specified args
	"""
	cfg = conf.RadialRepresentation
	gridType = cfg.function.__name__
	postfix = ["fd", "grid", gridType, "xmin%i" % cfg.xmin, 
			"xmax%i" % cfg.xmax, "count%i" % cfg.count, 
			"include_left%s" % cfg.include_left_boundary, 
			"include_right%s" % cfg.include_right_boundary]

	return postfix


@RegisterAll
def GetAngularPostfix(conf):
	"""
	Returns a "unique" list of strings string identifying the angular grid
	implied by the specified args
	"""
	cfg = conf.AngularRepresentation

	lmax = max([l1 for l1, l2, L, M in cfg.index_iterator])
	Llist = unique([L for l1, l2, L, M in cfg.index_iterator])
	Mlist = unique([M for l1, l2, L, M in cfg.index_iterator])

	def getSortedListString(l):
		if len(l) == 1:
			string = "%i" % l[0]
		else:
			if (diff(l) == 1).all():
				string = "%i-%i" % (l[0], l[-1]+1)
			else:
				string = "%s" % ("_".join(map(str, l)))

		return string

	postfix = ["angular"]
	postfix += ["lmax%i" % lmax]
	postfix += ["L%s" % getSortedListString(Llist)]	
	postfix += ["M%s" % getSortedListString(Mlist)]	

	return postfix


def GetPropagationPostfix(conf):
	"""
	Returns a "unique" list of strings identifying the propagation parameters
	implied by the specified PyProp config, 'conf'.
	"""
	cfg = conf.Propagation

	postfix = [
		"duration%s" % cfg.duration, 
		"dt%s" % cfg.timestep]

	return postfix

