"""
Description here
"""
__all__ = ["SingleParticleStates"]

import helium.namecontroller.namegenerator.AnalysisNameGenerator as NameGen
from __future__ import with_statement

class States:
	pass

class SingleParticleStates(States):
	"""
	Purpose: provide single particle states data

	get names, load data, nice interface
	"""

	def __init__(self, model, stateFilter, conf, ordering="ln"):
		self.Config = conf
		self.Model = model
		self.StateFilter = stateFilter
		self.Ordering = "ln"
		self.__States =  {}
		self.__Energies = {}

		#Generate single particle states file name
		self.FileName = NameGen.GetSingleParticleStatesFilenames(model=model, **args)

		#Load single particle states
		self.__Load(ordering)


	def GetEnergyFilteredRadialStates(self, l, energyFilter):
		"""
		Return all radial states corresponding to angular momentum 'l'
		and matching energy conditions in 'energyFilter'.
		"""
		getStates = lambda E, V: array([ V[i,:] for i in range(V.shape[0]) if energyFilter(E[i]) ])
		getEnergies = lambda E: filter(energyFilter, E)
		filteredStates = map(getStates, singleEnergies, singleStates)
		filteredEnergies = map(getEnergies, singleEnergies)

		return filteredEnergies, filteredStates	

	
	def IterateRadialStates(self, filter=lambda x: True):
		raise NotImplementedError("Not implemented yet!")

	
	def __Load(singleStatesFile):
		"""
		Load single particle states from file.
		"""
		def loadStates(f):
			#Create BSpline object
			cfgSection = pyprop.Section("RadialRepresentation", f.root.RadialEig._v_attrs.configObject)
			bspl = pyprop.BSPLINE()
			bspl.ApplyConfigSection(cfgSection)

			#Setup grid to check for correct phase convention
			phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
			phaseBuffer = zeros(2, dtype=complex)

			try:
				lList = f.root.RadialEig.l[:]
				for l in lList:
					node = f.getNode("/RadialEig/L%03i" % l)
					E = node.eigenvalues[:]
					V = node.eigenvectors[:]

					#assure correct phase convention (first oscillation should start out real positive)
					for i, curE in enumerate(E):
						bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
						phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
						V[:,i] *= exp(-1.0j * phase)
					
					#Store states and energies on object
					if ordering == "ln":
						self.__States[l] = transpose(V)
						self.__Energies[l] = E
					#elif ordering == "nl":
					#	self.__States[i] = zeros((), dtype='complex')
					#	for i, curE in enumerate(E):
					#		self.__States[i][l,:] = V[:,i]
					#		self.__Energies[i][l] = E[i]
						else:
							raise Exception("Unknown ordering: %s" % ordering)

		with f as tables.openFile(self.FileName, "r"):
			loadStates(f)

