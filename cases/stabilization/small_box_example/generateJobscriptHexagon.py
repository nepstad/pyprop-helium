def generateJobscript(pyTask,**args):
	"""
	generateJobscript(pyTask,**args)

	Generates a jobscript for stallo.uit.no 

	Parametres
	----------
	pyTask			: python script (and function) to be executed (string).
					Ex: '-c "execfile(\'small_box.py\'); ManyPropagations()'

	*TotalNumbProc*	: total number of processors, default 1
	*time*			: walltime, default [1,0,0] ([h,min,sec])
	*memory*		: memory, default 1000 (MB)
	*ProcPerNode*	: number of processors per node.
	*filename'		: name of jobscript, default autogenerated from processor info.
	
	
	"""

	# Userspesific parameters
	userName = 'sso059'
	projectAccount = 'fysisk' 



	if 'TotalNumbProc' in args:
		TotalNumbProc = args['TotalNumbProc']
	else : TotalNumbProc = 1

	if 'time' in args:
		time = args['time']
	else : time = [1,0,0]

	if 'memory' in args:
		memory = args['memory']
	else : memory = 1000

	if 'ProcPerNode' in args:
		ProcPerNode = args['ProcPerNode']
	else : ProcPerNode = 4

	if 'pyFunction' in args:
		pyFunction = args['pyFunction']
	else :
		pyFunction = ''

	if 'filename' in args:
		filename = args['filename']
	else :
		filename = 'job_' + str(lnodes) + '_' + str(ppn) + '_' + str(TotalNumbProc) + '_'	



	fileHandle = open(filename,'w')
	
	fileHandle.write('#!/bin/sh - \n' )
	fileHandle.write('#PBS -S /bin/sh \n')
	fileHandle.write('#PBS -N \"' + filename + '\"\n')
	
	
	wallTimeString = '#PBS -l walltime=' + str('%02.d' % time[0]) + ':'\
		+ str('%02.d' % time[1]) + ':' + str('%02.d' % time[2]) + ',mppwidth=' + str(TotalNumbProc) + str('\n')
	fileHandle.write(wallTimeString)	
	
	memString = '#PBS -l mppnppn=' +str(ProcPerNode) +  ',mppmem=' + str(memory) + 'mb \n'
	fileHandle.write(memString)
	
	fileHandle.write('#PBS -A ' + projectAccount + '\n')

	fileHandle.write('cd /work/' + userName + '/pyprop-helium/cases/stabilization/proper_calculations \n')
	
	mpirunString = 'aprun -n ' + str(TotalNumbProc) + ' -N ' + str(ProcPerNode) + ' python ' + pyTask + '\n'
	fileHandle.write(mpirunString)	
	
	print "The file " + filename + " has been generated."


global JobCount
JobCount = 0

def setJobCounter(startFrom):
	"""
	(Re)set JobCounter to a chosen start value.
	"""
	global JobCount
	JobCount = startFrom



def generateJobscriptWithPropagationParams(**args):
	"""
	generateJobscriptWithPropagationParams(**args)

	Generates a jobscript which executes the function ManyPropagations() in 
	small_box.py with the preferred parameters. 

	Parametres
	----------
	*freq* 	: frequency (string), defualt 'r_[10]'
	*E_0s* 	: electric field strength (string), default 'r_[30]'
	*cyc*	: cycles (string), default '[6]'
	*xsize*	: xsize (string), default '30'
	*xmax*	: xmax (string), default '30'
	*lmax* 	: lmax (string), default '5'
	*Ls*	: Ls (string), default '[0,1,2,3]'
	
	Notes
	-----
	Input parameters are strings.
	The function calls generateJobscript() and may pass arguments to this function.

	"""

	if 'freq' in args:
		freq = args['freq']
	else :
		freq='r_[10]'

	if 'E_0s' in args:
		E_0s = args['E_0s']
	else :
		E_0s = 'r_[30]'

	if 'cyc' in args:
		cyc = args['cyc']
	else :
		cyc = '[6]'
	
	if 'xsize' in args:
		xsize = args['xsize']
	else :
		xsize = '30'

	if 'xmax' in args:
		xmax = args['xmax']
	else :
		xmax = '30'

	if 'lmax' in args:
		lmax = args['lmax']
	else :
		lmax = '5'

	if 'Ls' in args:
		Ls = args['Ls']
	else :
		Ls = '[0,1,2,3]'

	global JobCount


	argString = freq + ',' + E_0s + ',' + cyc + ',' + xsize + ',' + xmax + ',' + Ls + ',' + lmax

	fileName = 'job_prop_' + str(JobCount) + '_'
	task = '-c "execfile(\'small_box.py\'); ManyPropagations(' + argString + ')"'
	
	generateJobscript(pyTask=task, filename=fileName,**args)

	#Increment counter
	JobCount += 1



def generateJobscriptWithGridParams(**args):
	"""
	generateJobscriptWithGridParams(**args)

	Generates a jobscript which executes the function FindBoundstatesAfterCheck() in 
	SandSfunctions.py with the preferred parameters. 

	Parametres
	----------
	*xsize*	: xsize (string), default '30'
	*xmax*	: xmax (string), default '30'
	*order*	: B-spline order
	*lmax* 	: lmax (string), default '5'
	*Ls*	: Ls (string), default '[0,1,2,3]'
	
	Notes
	-----
	Input parameters are strings.
	The function calls generateJobscript() and may pass arguments to this function.

	"""

	
	if 'xsize' in args:
		xsize = args['xsize']
	else :
		xsize = '30'

	if 'xmax' in args:
		xmax = args['xmax']
	else :
		xmax = '30'

	if 'order' in args:
		order = args['order']
	else :
		order = 5


	if 'lmax' in args:
		lmax = args['lmax']
	else :
		lmax = '5'

	if 'Ls' in args:
		Ls = args['Ls']
	else :
		Ls = '[0,1,2,3]'


	global JobCount


	argString = xsize + ',' + xmax + ',' + order +',' + Ls + ',' + lmax

	fileName = 'job_bound_' + str(JobCount) + '_'
	task = '-c "execfile(\'SandSfunctions.py\'); FindBoundstatesAfterCheck(' + argString + ')"'
	
	generateJobscript(pyTask=task, filename=fileName,**args)

	#Increment counter
	JobCount += 1



