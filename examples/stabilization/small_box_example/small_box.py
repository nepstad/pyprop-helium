from __future__ import with_statement
import sys
import logging

sys.path.append("..")
import stabilization

import helium.propagation.tasks as Tasks
from helium.eigenvalues.eigenvalues import FindEigenvaluesInverseIterationsPiram
from helium.eigenvalues.eigenvalues import SaveEigenvalueSolverShiftInvert
from helium.propagation.propagate import Propagate
from helium.eigenvalues.eigenstates import LoadBoundstateIndex
from helium.analysis.observables import ContinuumObservables
import pyprop
from numpy import *

#set logging level
logging.basicConfig(level=logging.DEBUG)

def FindBoundstates(lmax,Ls,xsize,xmax,order):
	"""
	Find some bound states of Helium using Piram + inverse iterations.
	The eigenpairs are also stored in a file (name is autogenerated).
	
	Parameters
	----------
	Ls : list of large Ls.
	
	Returns a list eigenvalues.
	"""

	def UpdateConfig(conf):
		
		return pyprop.Config(conf.cfgObj)


	#Load a config file
	conf = pyprop.Load("eigenvalues.ini")

	#Change config according to parametres
	conf.SetValue("RadialRepresentation", "xsize", xsize)
	conf.SetValue("RadialRepresentation", "xmax", xmax)
	conf.SetValue("RadialRepresentation", "order", order)
	
	
	conf = UpdateConfig(conf)


	#Setup solver and run eigenvalue iterations for L = 0,1
	eigs = {}
	for L in Ls:
		logging.info("Calculating bound states for L = %i..." % L)
		angIt = conf.AngularRepresentation.index_iterator
		angIt.L = [L]
		angIt.lmax = lmax
		conf.SetValue("AngularRepresentation", "index_iterator", angIt)
		with pyprop.EnableRedirect():
			solver, gmresSolver, E = FindEigenvaluesInverseIterationsPiram(conf)

		#Save eigenpairs to HDF-file (name is autogenerated)
		SaveEigenvalueSolverShiftInvert(solver, gmresSolver)
		
		eigs['L_%i' % L] = E

	return eigs


def SetupTasks():
	"""
	Setup some propagation tasks defined in the Helium module.

	Returns a list of propagation tasks.
	"""
	tasks = []
	tasks += [Tasks.SaveWavefunction(True)]
	tasks += [Tasks.ProgressReport()]
	#tasks += [Tasks.DisplayGMRESError()]
	
	return tasks


def PropagationExecutor(frequency = 5.0, E_0 = 1.0, cycles = 6, xsize = 60, xmax = 30.0, Ls= 3, lmax = 5, **args):
	"""
	Calculate the effect of a short, intense laser pulse on the Helium
	ground state, by time propagation using the Cayley propagator.

	Uses the 'propagation.ini' config file.

	Returns helium.propagate.Propagate object.

	"""

	if propagationConfig in args:
		propagationConfig = args['propagationConfig']
	else : 
		propagationConfig = 'propagation.ini'	

	def UpdateConfig(conf):
		
		return pyprop.Config(conf.cfgObj)

	#Load config
	conf = pyprop.Load(propagationConfig)

	#Change config according to parametres
	conf.SetValue("Names", "output_file_name", "data/omega_%.2f_E_%.2f_cyc_%i_size_%i/proper_propagation.h5"%(frequency,E_0,cycles,xsize))
	conf.SetValue("PulseParameters", "frequency", frequency)
	conf.SetValue("PulseParameters", "amplitude", E_0/frequency)
	conf.SetValue("PulseParameters", "cycles", cycles)
	conf.SetValue("RadialRepresentation", "xsize", xsize)
	conf.SetValue("RadialRepresentation", "xmax", xmax)
	
	angIt = conf.AngularRepresentation.index_iterator
	angIt.L = Ls
	angIt.lmax = lmax
	conf.SetValue("AngularRepresentation", "index_iterator", angIt)




	conf = UpdateConfig(conf)
	
	#Setup propagation tasks
	tasks = SetupTasks()
	#Setup propagate object
	prop = Propagate(conf, tasks, 50)
	
	#Run propagation
	prop.run()
	
	return prop

def test(frequencies, E_0s, cycles, xsize,xmax, Ls, lmax):
		print "HELLO TROMSO!"
		print frequencies
		print E_0s[0]
		print cycles, xsize, xmax
		print Ls[1], lmax

def ManyPropagations(frequencies, E_0s, cycles, xsize,xmax, Ls, lmax,**args):
	"""
	Loops over PropegationSmallBox() for the freq. and E_0 given as parametres.

	Parametres
	----------
	frequencies : Array, containing the frequencies for which we want to propagate the problem.
	E_0s : Array, containing the electric field strengths for which we want to propagate the problem.
	"""
	print "HELLO TROMSO"	

	config_list = []
	for freq in frequencies:
		for E_0 in E_0s:
			for cyc in cycles:
				prop = PropagationExecutor(frequency = freq, E_0 = E_0, cycles = cyc, 
						xsize = xsize, xmax = xmax, Ls= Ls, lmax = lmax,**args)
				del prop


def CalculateIonization(conf):
	"""
	Calculates total-, single- and double ionization probabilities for
	a propagated wavefunction, defined by the Pyprop config object 'conf'.
	
	Returns a 'DoubleContinuumObservables' instance.
	"""

	#Setup/load wavefunction
	logging.info("Load wavefunction...")
	with pyprop.EnableRedirect():
		psiFileName = conf.Names.output_file_name
		psi = pyprop.CreateWavefunctionFromFile(psiFileName)

	#Observables
	logging.info("Calculating observables...")
	#isIonized = lambda E: E > 0.0
	observables = DoubleContinuumObservables(psi, conf)
	observables.Setup()
	logging.info("  Calculation double ionization...")
	
	doubleIon = observables.GetDoubleIonizationProbability()
	totalIon = observables.TotalIonizationProbability
	logging.info("    Total ionization = %s" % totalIon)
	logging.info("    Double ionization = %s" % doubleIon)
	
	return observables

def CalculateTotalIonization(conf):
	"""
	Calculates total ionization probabilities for
	a propagated wavefunction, defined by the Pyprop config object 'conf'.
	
	Returns a 'DoubleContinuumObservables' instance.
	"""

	#Setup/load wavefunction
	logging.info("Load wavefunction...")
	with pyprop.EnableRedirect():
		psiFileName = conf.Names.output_file_name
		psi = pyprop.CreateWavefunctionFromFile(psiFileName)

	#Observables
	logging.info("Calculating observables...")
	#isIonized = lambda E: E > 0.0
	observables = ContinuumObservables(psi, conf)
	observables.Setup()
	logging.info("  Calculation double ionization...")
	
	totalIon = observables.IonizationProbability
	logging.info("    Total ionization = %s" % totalIon)
	
	return observables
