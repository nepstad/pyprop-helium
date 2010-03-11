"""
NameGenerator
=============

Generate names specific to stabilization cases

"""

def GetFinalOutputFilename(conf):
	customPrefix = "wavefunction"
	params = conf.PulseParameters

	#get salient parameters
	freq = str(params.frequency)
	duration = "%3.1f" % params.duration
	timestep = "%1.e" % conf.Propagation.timestep
	pulseName = params.time_function.__name__
	amplitude = str(params.amplitude)
	
	#construct string representation of parameters
	filename = "_".join([customPrefix, pulseName, "I", amplitude, "freq", freq, "dt", timestep, "duration", duration])
	
	return filename

