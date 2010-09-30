import os, sys
import commands
from datetime import timedelta
sys.path.append("..")
sys.path.append("../pyprop")
from helium.configtools import UpdateConfig

INSTALLATION = getattr(os.environ, "INSTALLATION", "N/A")
if INSTALLATION == "hexagon":
	from pyprop.utilities import submitpbs_hexagon
elif INSTALLATION == "stallo":
	from pyprop.utilities import submitpbs_stallo
	
#------------------------------------------------------------------------------------
#                       Job Submit Functions
#------------------------------------------------------------------------------------
def RunSubmitFullProcCount(function, account, executable, conf, updateParams, **args):
	angCount = GetAngularBasisSize(conf, updateParams)
	return RunSubmit(function, account, executable, conf, angCount, procPerNode, updateParams, **args)


def RunSubmit(function, account, executable, conf, procCount=1, procPerNode=4, procMemory="1000M", \
			walltime=timedelta(minutes=30), depend=None, writeScript=False, scriptname="test.job", \
			**args):
	"""
	Runs a function on the compute nodes.
	
	function is either a string (name of the function) or a function declaration. 

	function is run by creating a job script for run-function.py, and passing the
	function name and all arguments to the job script.
	"""

	if isinstance(function, str):
		arg1 = function
	else:
		arg1 = function.func_name

	jobId = None
	if INSTALLATION == "hexagon":
		submit = submitpbs_hexagon.SubmitScript()
		submit.procs = procCount
		submit.ppn = min(procPerNode, procCount)
		submit.proc_memory = procMemory
		submit.executable = executable
		submit.parameters = arg1 + 
		submit.walltime = walltime
		submit.account = account
		submit.depend = depend
		submit.WriteScript(scriptname)
		if not writeScript:
			jobId = submit.Submit()
		else:
			print "Ok, just writing the submit script then"
			
	elif INSTALLATION == "stallo":
		submit = submitpbs_stallo.SubmitScript()
		submit.procs = procCount
		submit.ppn = min(procPerNode, procCount)
		submit.proc_memory = procMemory
		submit.executable = executable
		submit.parameters = updateParams
		submit.walltime = walltime
		submit.interconnect = "ib"
		submit.WriteScript("test.job")
		if not writeScript:
			jobId = submit.Submit()
		else:
			print "Ok, just writing the submit script then"

	elif INSTALLATION == "local":
		raise Exception("please to implement")
	
	else:
		raise Exception("Unknown installation '%s'" % INSTALLATION)

	return jobId

#------------------------------------------------------------------------------------
#                       Job Submit Functions
#------------------------------------------------------------------------------------
def GetAngularBasisSize(conf, params):
	"""
	Returns the number of spherical harmonic basis functions for a 
	given argument list
	"""
	conf = UpdateConfig(conf, params)
	return len([1 for i in conf.AngularRepresentation.index_iterator])

