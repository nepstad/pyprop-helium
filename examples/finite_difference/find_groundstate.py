from numpy import diff
import sys
sys.path.append("../../pyprop")
sys.path.append("../..")
import pyprop

import helium
import helium.core
import helium.core.preconditioner
from helium.propagation.propagate import Propagate
from helium.propagation.tasks import ProgressReport, SaveWavefunction
from helium.propagation.tasks import PropagationTask
from helium.utils import UpdatePypropProjectNamespace
UpdatePypropProjectNamespace(pyprop.ProjectNamespace)


class EnergyReport(PropagationTask):
	def __init__(self):
		pass

	def setupTask(self, prop):
		pass

	def preprocess(self, prop):
		pass

	def callback(self, prop):
		curEnergy =  prop.GetEnergyExpectationValue()
		pyprop.PrintOut("Energy = %s" % curEnergy)

	def postProcess(self, prop):
		pass


def FindGroundstate():
	"""Find groundstate using imaginary time propagation
	
	"""
	timers = pyprop.Timers()

	timers["Setup"].Start()

	#Load config
	conf = pyprop.Load("groundstate.ini")

	#Initialize propagation tasks
	tasks = [ProgressReport(), EnergyReport(), SaveWavefunction(False)]

	#Setup propagation
	prop = Propagate(conf, tasks, 20)

	timers["Setup"].Stop()

	print "psi shape = ", prop.Problem.psi.GetData().shape
	dx = diff(prop.Problem.psi.GetRepresentation().GetGlobalGrid(0))[0]
	pyprop.PrintOut("Radial grid spacing = %1.16f" % dx)

	#Run propagation
	timers["Propagate"].Start()
	prop.run()
	timers["Propagate"].Stop()

	if pyprop.IsMaster():
		print
		print timers

	return prop


if __name__ == "__main__":
	_ = FindGroundstate()
