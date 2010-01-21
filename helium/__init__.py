"""
Helium
======

The Helium package provides fundamental functionality related to Helium-type problems,
such as potentials, analysis and eigenpair calculations.
"""
__all__ = ["analysis", "core", "eigenvalues", "namecontroller", "siteconfig", "utils"]

import sys
import siteconfig as SiteConfig

#Add Pyprop location to path
sys.path.append(SiteConfig.PypropLocation)

#Contains function and class name from package that might appear in config files
#To be passed to Pyprop for resolving during loading of said files.
from utils import ProjectNamespace

#import pyprop
#import namecontroller as NameController
#import analysis as Analysis
#import utils as Utils

