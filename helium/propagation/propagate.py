"""
Propagate
=========

Classes which wraps Pyprop problems and structures the propagation flow into
four distinct steps: init, pre-process, run and post-process. This has turned
out to be a convenient way to treat most propagation problems, as most operations
one might wish to perform invariably fall into one of said categories.

The general idea here is that all operations besides the actual propagation
are implemented as propagation tasks, which are initialized by the user and
then passed as a list to the Propagate class (or derivatives thereof) - see
tasks.py for examples of how these should/could be implemented.  

"""

import pyprop
from helium.utils import GetClassLogger

class Propagate(object):
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
		for _ in self.Problem.Advance(self.NumberOfCallbacks):
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
	
	
class PropagateYieldCondition(Propagate):
	"""
	Setup and run a Pyprop problem, with a given condition of when to
	return flow control to user. The condition is specified in terms 
	of a function which takes the current state of the system as an
	argument (a Pyprop problem).
	
	See also Propagate.
	
	"""
	
	def __init__(self, conf, propagationTasks, yieldCondition):
		#Call on constructor of superclass
		super(PropagateYieldCondition, self).__init__(conf, propagationTasks, None)
		
		#Keep the yield condition (function)
		self.YieldCondition = yieldCondition
		
	def run(self):
		"""
		Propagate problem until end time. Check if YieldCondition is fulfilled
		in each time step, if so, call all propagation tasks.
		"""
		for _ in self.Problem.Advance(True):
			if self.YieldCondition(self.Problem):
				for task in self.PropagationTasks:
					task.callback(self.Problem)
		
		#run postprocessing
		self.postProcess()
		
		

		