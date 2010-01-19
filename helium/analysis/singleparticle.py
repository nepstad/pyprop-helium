"""
SingleParticle
==============

Basic single particle states operations are provided. These are:

	1)Generate file names
	2)Load single particle data
	3)?

Provides:

"""
from __future__ import with_statement

__all__ = ["SingleParticleStates"]

import helium.namecontroller.namegenerator as NameGen
import tables
import pyprop
from numpy import array, zeros, arctan2, exp, transpose, double, complex, imag, real

class States:
	pass

class SingleParticleStates(States):
	"""
	Purpose: provide single particle states data

	get names, load data, nice interface
	"""

	def __init__(self, model, energyFilter, conf, ordering="ln"):
		"""
		"""
		self.Config = conf
		self.Model = model
		self.EnergyFilter = energyFilter
		self.Ordering = "ln"
		self.__States =  {}
		self.__Energies = {}

		#Generate single particle states file name
		self.FileName = NameGen.GetSingleParticleStatesFilename(conf, model)

		#Load single particle states
		self.__Load()


	def GetEnergyFilteredRadialStates(self, l):
		"""
		Return all radial states corresponding to angular momentum 'l'
		and matching energy conditions in self.StateFilter.

		Input
		-----
		l: angular momentum

		Returns: energies (1D array) and eigenstates (2D array)
		"""
		numStates = self.__States[l].shape[1]
		E = self.__Energies[l]
		V = self.__States[l]
		filteredStates = array([V[:,i] for i in range(numStates) if self.EnergyFilter(E[i])])

		filteredEnergies = filter(self.EnergyFilter, self.__Energies[l])

		return filteredEnergies, filteredStates.transpose()

	
	def IterateRadialStates(self):
		"""
		"""
		for l in self.__States.iterkeys():
			E,V = self.GetEnergyFilteredRadialStates(l)
			yield int(l), V


	def __Load(self):
		"""
		Load single particle states from file.
		"""
		def loadStates(f):
			#Create BSpline object
			cfgSection = self.Config.RadialRepresentation
			bspl = pyprop.BSPLINE()
			bspl.ApplyConfigSection(cfgSection)

			#Setup grid to check for correct phase convention
			phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
			phaseBuffer = zeros(2, dtype=complex)

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
				if self.Ordering == "ln":
					self.__States[l] = V
					self.__Energies[l] = E
				#elif ordering == "nl":
				#	self.__States[i] = zeros((), dtype='complex')
				#	for i, curE in enumerate(E):
				#		self.__States[i][l,:] = V[:,i]
				#		self.__Energies[i][l] = E[i]
				else:
					raise Exception("Unknown ordering: %s" % self.Ordering)

		#Do the loading; always close file
		with tables.openFile(self.FileName, "r") as f:
			loadStates(f)

