"""
Helium
======

The Helium package provides fundamental functionality related to Helium-
type problems, such as potentials, analysis and eigenpair calculations.

"""
import sys
import os

#Add Pyprop location to path
PypropLocation = os.path.realpath("%s/../pyprop" % \
		os.path.split(os.path.realpath(__file__))[0])
sys.path.append(PypropLocation)

#Contains function and class name from package that might appear in config
#files. To be passed to Pyprop for resolving during loading of said files.
from utils import ProjectNamespace

__all__ = ["analysis", "core", "eigenvalues", "namecontroller", "siteconfig", \
		"utils", "configtools"]
