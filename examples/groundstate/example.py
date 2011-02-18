#!/bin/env python
import sys
sys.path.append("../..")

import helium
import helium.core
import helium.core.preconditioner
from helium.eigenvalues.eigenvalues import FindEigenvaluesInverseIterationsPiram
import pyprop

#Update Pyprop project namespace with objects from Helium project namespace.
#This is required in order for Pyprop to resolve Helium classes when loading
#config files.
from helium.utils import UpdatePypropProjectNamespace
UpdatePypropProjectNamespace(pyprop.ProjectNamespace)


def FindBoundstatesExample():
	"""Find the 5 lowest eigenvalues of Helium using shift-invert IRAM

	"""
	#Load config 
	conf = pyprop.Load("eigenvalues.ini")

	#Compute eigenvalues
	solver, siSOlver, eigenVals = FindEigenvaluesInverseIterationsPiram(conf)

	#Print eigenvalues
	pyprop.PrintOut("\nFound eigenvalues:")
	pyprop.PrintOut(" ".join(["%f" % E for E in eigenVals.real]))


if __name__ == "__main__":
	FindBoundstatesExample()

