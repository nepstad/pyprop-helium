"""
eigenstates
===========

Basic eigenstate operations are provided. These are:

	1)Load/store eigenstates from disk
	2)Load (single-L) eigenstate into multi-L wavefunction


"""
from __future__ import with_statement
import tables
from numpy import r_, s_, array, nonzero, ravel
import pyprop
from helium.utils import RegisterAll, RegisterProjectNamespace
from helium.namecontroller import namegenerator as NameGen
from helium.configtools import GetL
from helium.namecontroller.namegenerator import GetEigenvectorDatasetName
from helium.namecontroller.namegenerator import GetEigenvectorDatasetPath
from helium.utils import GetClassLogger, GetFunctionLogger

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

		return boundEnergies


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


@RegisterProjectNamespace
@RegisterAll
def LoadBoundstateIndex(psi, filename, eigenvectorIndex):
	"""
	Loads an eigenstate from a file <filename> into a wavefunction <psi>. 
	This function works in parallel.

	Because eigenstates can be found independently for each "L", 
	the number and order of angular indices are different in the
	stored eigenstates, from that of the wavefunction to load the 
	eigenstate into.

	The eigenstate config file is used to create an array of CoupledIndex elements
	to which the local CoupledIndex elements from the wavefunction are compared

	"""
	angularRank = 0
	logger = GetFunctionLogger()

	logger.info("Loading bound state #%i from %s" % (eigenvectorIndex, filename))
	
	with tables.openFile(filename, "r") as f:
		#Get eigenvector dataset
		eigVec = f.getNode(GetEigenvectorDatasetPath(eigenvectorIndex))

		#Get corresponding eigenvalue
		E = f.root.Eig.Eigenvalues[eigenvectorIndex]

		#find angular indices of the saved eigenstate
		conf = pyprop.LoadConfigFromFile(filename, "/Eig")
		eigCoupledIndices = array([i for i in conf.AngularRepresentation.index_iterator])

		#find local angular indices of the wavefunction to load into
		repr = psi.GetRepresentation()
		angRepr = repr.GetRepresentation(angularRank)
		localIndices = array(repr.GetLocalGrid(angularRank), dtype=int)

		#Load the local angular indices into the wavefunction
		psi.Clear()
		for curLocalIdx, curGlobalIdx in enumerate(localIndices):
			#Get the CoupledIndex to this angular index
			curCoupledIndex =  angRepr.Range.GetCoupledIndex(int(curGlobalIdx))

			#Find out the index in the eigenvector corresponding
			eigIdx, = nonzero(ravel(eigCoupledIndices == curCoupledIndex))
			if len(eigIdx) == 0:
				continue
			elif len(eigIdx) > 1:
				raise Exception("What? %s" % eigIdx)
			else:
				eigIdx = eigIdx[0]
				

			#Slice the wavefunction at the correct angular index
			psiSlice = [s_[:]]*psi.GetRank()
			psiSlice[angularRank] = curLocalIdx

			#Slice the eigenstate at the correct angular index
			eigSlice = [s_[:]]*psi.GetRank()
			eigSlice[angularRank] = eigIdx

			psi.GetData()[curLocalIdx,:,:] = eigVec[eigIdx,:,:]

	return E

