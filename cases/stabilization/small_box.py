execfile("init.py")

def SetupTasks():
	tasks = []
	tasks += [helium.propagation.tasks.SaveWavefunction(True)]
	tasks += [helium.propagation.tasks.ProgressReport()]
	tasks += [helium.propagation.tasks.DisplayGMRESError()]
	
	return tasks


def TestPropagationSmallBox():
	conf = pyprop.Load("propagation.ini")
	tasks = SetupTasks()
	prop = helium.propagation.propagate.Propagate(conf, tasks, 50)
	prop.run()

	return prop
