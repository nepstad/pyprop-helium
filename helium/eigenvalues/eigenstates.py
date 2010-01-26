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
import pyprop
from ..utils import RegisterAll
from ..namecontroller import namegenerator as NameGen

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
		E = self._Energies[L]
		boundIdx = filter(lambda i: E[i] < ionizationThreshold, r_[:len(E)])

	
	def _Load(self):
		for filename in self.FileNames:
			L = 0
			eigPsi = pyprop.CreateWavefunctionFromFile(filename, GetEigenvectorDatasetPath(0))

			curPsiList = []
			curEnergyList = []

			#Get energies
			with tables.openFile(filename, "r") as f:
				eigenvalues = f.root.Eig.Eigenvalues[:]
			
			for i in len(eigenvalues):
				psi = eigPsi.Copy()
				pyprop.serialization.LoadWavefunctionHDF(filename, GetEigenvectorDatasetPath(i), psi)
				curPsiList.append(psi)
				curEnergyList.append(eigenvalues[i])

			self._States[L] = curPsiList
			self._Energies[L] = curEnergyList
