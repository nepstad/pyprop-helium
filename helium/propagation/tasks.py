"""
PropagationTasks
================

Defines standard tasks to be performed during propagation

"""
from __future__ import with_statement
import os.path
import time
import tables
from pyprop import PrintOut
from helium.utils import RegisterAll, GetClassLogger


def CreatePath(absFileName):
	filePath = os.path.dirname(absFileName)
	if not os.path.exists(filePath):
		os.makedirs(filePath)
	

class PropagationTask:
	def __init__(self):
		pass
	
	def setupTask(self, prop):
		pass

	def callback(self, prop):
		pass

	def postProcess(self, prop):
		pass


@RegisterAll
class ProgressReport(PropagationTask):
	"""
	Print some progress information during propagation, and store it
	in a HDF5-file when finished (time, norm and projection on initial state)
	"""

	def __init__(self):
		self.Logger = GetClassLogger(self)
		self.TimeList = []
		self.NormList = []
		self.CorrList = []
		self.StartTime = -1
		self.InitialPsi = None

	def setupTask(self, prop):
		self.Logger.info("Setting up task...")
		self.StartTime = time.time()
		self.InitialPsi = prop.psi.Copy()
		self.OutputFileName = prop.Config.Names.output_file_name

		#check if output dir exist, create if not
		CreatePath(self.OutputFileName)

	def callback(self, prop):
		t = prop.PropagatedTime
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(self.InitialPsi))**2
		eta = self.__EstimateETA(prop)
		self.TimeList += [t]
		self.NormList += [norm]
		self.CorrList += [corr]
		FormatDuration = lambda t: time.strftime("%Hh %Mm %Ss", time.gmtime(t))
		PrintOut("t = %.2f / %.2f; N = %.15f; Corr = %.12f, ETA = %s" % (t, prop.Duration, norm, corr, FormatDuration(eta)))

	def postProcess(self, prop):
		"""
		Store problem information collected during propagation
		"""
		if pyprop.ProcId == 0:
			nodeNames = ["SampleTimes", "Norm", "InitialCorrelation"]
			with tables.openFile(self.OutputFileName, "r+") as h5file:
				for nodeName in nodeNames:
					if nodeName in h5file.root:
						f.removeNode(f.root, nodeName, recursive=True)

				h5file.createArray("/", "SampleTimes", self.TimeList)
				h5file.createArray("/", "Norm", self.NormList)
				h5file.createArray("/", "InitialCorrelation", self.CorrList)


	def __EstimateETA(self, prop):
		"""
		Estimates remaining time before propagation is finished
		"""
		curTime = time.time() - self.StartTime
		totalTime = (curTime / prop.PropagatedTime) * prop.Duration
		eta = totalTime - curTime
		return eta


class DisplayGMRESError(PropagationTask):
	"""
	Print GMRES solver error at each callback
	"""
	def callback(self, prop):
		PrintOut(prop.Propagator.Solver.GetErrorEstimateList())


class SaveWavefunction(PropagationTask):
	"""
	Save wavefunction after propagation, and, if specified, for each
	callback during propagation
	"""

	def __init__(self, storeDuringPropagation):
		self.StoreDuringPropagation = storeDuringPropagation
		self.Counter = 0

	def setupTask(self, prop):
		#get output filename
		self.OutputFileName = prop.Config.Names.output_file_name

		#check if output dir exist, create if not
		CreatePath(self.OutputFileName)

		#store the initial wavefunction
		prop.SaveWavefunctionHDF(self.OutputFileName, "/initialWavefunction")

	def callback(self, prop):
		if self.StoreDuringPropagation:
			#create unique filename
			filename = "%s_%03i.h5" % (self.OutputFileName.strip(".h5"), self.Counter)
			
			#store current wavefunction and propagation time
			prop.SaveWavefunctionHDF(filename, "/wavefunction")
			with tables.openFile(filename, "r+") as h5:
				h5.setNodeAttr("/wavefunction", "prop_time", prop.PropagatedTime)

			self.Counter += 1

	def postProcess(self, prop):
		#store the final wavefunction
		prop.SaveWavefunctionHDF(self.OutputFileName, "/wavefunction")

