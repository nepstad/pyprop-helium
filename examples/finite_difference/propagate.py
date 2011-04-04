from numpy import sin, cos, pi
import sys
sys.path.append("../../pyprop")
sys.path.append("../..")
import pyprop

import helium
import helium.core
import helium.core.preconditioner
from helium.propagation.propagate import Propagate
from helium.propagation.tasks import ProgressReport, SaveWavefunction
from helium.utils import UpdatePypropProjectNamespace
UpdatePypropProjectNamespace(pyprop.ProjectNamespace)

import serialization
from serialization import LoadGroundstate


def RunPropagation(configFile="propagate.ini"):
	"""Perform propagation defined by 'config.ini' file.

	"""
	#Load config
	conf = pyprop.Load(configFile)

	#Setup propagation tasks
	tasks = [ProgressReport(), SaveWavefunction(False)]

	#Setup propagation
	prop = Propagate(conf, tasks, 100)

	#Run propagation
	prop.run()


	return prop


def LaserFunctionSimpleLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= cos(t * conf.frequency);
	else:
		curField = 0
	return curField
#Put laser function in pyprop project namespace so that config files are
#loaded properly.
pyprop.ProjectNamespace["LaserFunctionSimpleLength"] = LaserFunctionSimpleLength


if __name__ == "__main__":
	_ = RunPropagation()
