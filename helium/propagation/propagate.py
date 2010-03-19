"""
Propagate
=========

"""

import pyprop

class Propagate:
	"""
	Setup and run a Pyprop problem.
	
	Propagate wavefunction from	T_start to T_end. Also perform
	all given PropagationTasks during the propagation phase.

	"""
	def __init__(self, conf, propagationTasks, numberOfCallbacks):
		self.PropagationTasks = propagationTasks
		self.Config = conf
		self.NumberOfCallbacks = numberOfCallbacks

		#setup Pyprop problem from config
		self.Problem = pyprop.Problem(self.Config)
		self.Problem.SetupStep()

		#run pre-propagation step for all tasks
		for task in self.PropagationTasks:
			task.setupTask(self.Problem)

		en = self.Problem.GetEnergyExpectationValue()
		print "Initial state energy = %s" % en 
			

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
