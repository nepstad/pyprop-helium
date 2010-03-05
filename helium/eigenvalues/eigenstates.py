"""
eigenstates
===========

Basic eigenstate operations are provided. These are:

	1)Generate file names
	2)Load/store eigenstates from disk
	3)?


"""
from __future__ import with_statement
import tables
from numpy import r_
import pyprop
from helium.utils import RegisterAll
from helium.namecontroller import namegenerator as NameGen
from helium.configtools import GetL
from helium.namecontroller.namegenerator import GetEigenvectorDatasetName
from helium.namecontroller.namegenerator import GetEigenvectorDatasetPath
@RegisterAll
class Eigenstates(object):
	def __init__(self, conf):
		"""
		"""
		self.Config = conf
		self._States =  {}
		self._Energies = {}

		#Generate bound state file names
		self.FileNames = NameGen.GetAllBoundstateFilenames(conf)

		#Load single particle states
		self._Load()


	def GetBoundstates(self, L, ionizationThreshold):
		"""Return list of bound state wavefunctions
		"""
		E = self._Energies[L]
		boundIdx = filter(lambda i: E[i] < ionizationThreshold, r_[:len(E)])
		boundStates = [self._States[L][idx] for idx in boundIdx]

		return boundStates


	def GetBoundEnergies(self, L, ionizationThreshold):
		"""Return list of bound state energies
		"""
		E = self._Energies[L]
		isBound = lambda en: en < ionizationThreshold
		boundEnergies = [curE for curE in E if isBound(E)]

		return self._Energies[L][boundIdx]


	def IterateBoundstates(self, ionizationThreshold):
		"""Iterate over all bound states
		"""
		isBound = lambda en: en < ionizationThreshold
		for L in self._Energies.iterkeys():
			for idx, E in enumerate(self._Energies[L]):
				if isBound(E): 
					psi = self._States[L][idx]
					yield L, E, psi

	
	def _Load(self):
		for filename in self.FileNames:
			#Get total angular momentum for this file
			with pyprop.EnableRedirect():
				conf = pyprop.LoadConfigFromFile(filename, '/Eig')
			L = GetL(conf)

			#Create eigenstate wavefunction
			with pyprop.EnableRedirect():
				eigPsi = pyprop.CreateWavefunctionFromFile(filename, \
					GetEigenvectorDatasetPath(0))

			curPsiList = []
			curEnergyList = []

			#Get energies
			with tables.openFile(filename, "r") as f:
				eigenvalues = f.root.Eig.Eigenvalues[:]
		
			#Load all eigenstates
			for i in range(eigenvalues.size):
				with pyprop.EnableRedirect():
					psi = eigPsi.Copy()
					pyprop.serialization.LoadWavefunctionHDF(filename, \
						GetEigenvectorDatasetPath(i), psi)
				curPsiList.append(psi)
				curEnergyList.append(eigenvalues[i])

			self._States[L] = curPsiList
			self._Energies[L] = curEnergyList

