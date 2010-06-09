"""
Propagate
=========

"""

import pyprop
from helium.utils import GetClassLogger

class Propagate:
	"""
	Setup and run a Pyprop problem.
	
	Propagate wavefunction from	T_start to T_end. Also perform
	all given PropagationTasks during the propagation phase.

	"""
	def __init__(self, conf, propagationTasks, numberOfCallbacks):
		self.Logger = GetClassLogger(self)
		self.PropagationTasks = propagationTasks
		self.Config = conf
		self.NumberOfCallbacks = numberOfCallbacks

		#setup Pyprop problem from config
		self.Problem = pyprop.Problem(self.Config)
		self.Problem.SetupStep()

		#run pre-propagation step for all tasks
		for task in self.PropagationTasks:
			task.setupTask(self.Problem)

		#Calculate intial state energy
		tmpPsi = self.Problem.psi.Copy()
		en = self.GetEnergyExpectationValue(self.Problem.psi, tmpPsi).real
		self.Logger.info("Initial state energy = %s" % en) 
			

	def run(self):
		"""
		Propagate problem until end time.
		"""
		for t in self.Problem.Advance(self.NumberOfCallbacks):
			for task in self.PropagationTasks:
				task.callback(self.Problem)
		
		#run postprocessing
		self.postProcess()

	
	def postProcess(self):
		"""
		Run postprocessing for all propagation tasks
		"""
		for task in self.PropagationTasks:
			task.postProcess(self.Problem)


	def GetEnergyExpectationValue(self, psi, tmpPsi):
		"""
		Calculates the total energy of the problem by finding the expectation value 
		of the time-independent part of the Hamiltonian. Assumes that Tensor potentials
		and BasisPropagator are used.
		"""
		energy = 0.0
		
		#Iterate over all potentials, check for time dependence
		for pot in self.Problem.Propagator.BasePropagator.PotentialList:
			if not pot.IsTimeDependent:
				energy += pot.GetExpectationValue(psi, tmpPsi, 0, 0)

		return energy
	