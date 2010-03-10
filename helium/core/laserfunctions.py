"""
LaserFunctions
==============

Supply time dependence for time-dependent potentials, such as laser fields.
"""

from ..utils import RegisterAll, RegisterProjectNamespace
from numpy import sin, cos, pi


@RegisterProjectNamespace
@RegisterAll
def LaserFunctionVelocity(conf, t):
	phase = 0.0
	if hasattr(conf, "phase"):
		phase = conf.phase
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= - cos(t * conf.frequency + phase);
	else:
		curField = 0
	return curField


@RegisterProjectNamespace
@RegisterAll
def LaserFunctionFlatTop(conf, t):
	pulseStart = 0
	if conf.Exists("pulse_start"):
		pulseStart = conf.pulse_start
	curField = conf.amplitude * cos(t * conf.frequency);

	if (t > conf.pulse_duration) or (t < pulseStart):
		curField = 0
	elif 0 <= t < conf.ramp_on_time:
		curField *= sin(t * pi / (2*conf.ramp_on_time))**2;
	elif t > conf.pulse_duration - conf.ramp_off_time:
		curField *= sin((conf.pulse_duration - t) * pi / (2*conf.ramp_off_time))**2;
	else:
		curField *= 1

	return curField


@RegisterProjectNamespace
@RegisterAll
def LaserFunctionSimpleLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= cos(t * conf.frequency);
	else:
		curField = 0
	return curField


@RegisterProjectNamespace
@RegisterAll
def LaserFunctionLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		T = conf.pulse_duration
		w = conf.frequency
		curField *= sin(pi*t/T)*(-2*pi/T * cos(pi*t/T) * cos(w*t) + w*sin(pi*t/T)*sin(w*t))
	else:
		curField = 0
	return curField


