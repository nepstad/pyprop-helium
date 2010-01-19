"""
Helium
======

Provides:

"""
__all__ = ["analysis", "core", "eigenvalues", "namecontroller", "siteconfig"]

import sys
import siteconfig as SiteConfig

#Add Pyprop location to path
sys.path.append(SiteConfig.PypropLocation)

import pyprop
import namecontroller as NameController


