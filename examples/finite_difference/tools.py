"""
Useful finite-difference related tools that do not require Pyprop

"""
from numpy import unique, array
import pylab

def FindProcNumber(N, maxProcPerRank=40):
	"""
	Determine possible Nproc = q x p for given number of grid points N

	Should have N % q == 0 and N % q == 0

	"""

	possibleProcs = [p*q for p in range(1,maxProcPerRank) for q in \
		range(1,maxProcPerRank) if N%p == 0 and N%q == 0] 

	possibleProcs = unique(possibleProcs)

	return possibleProcs



def ParseIPMOutput(filename):
	mpiTime = []
	mpiTotTime = []
	wallTime = []
	totTime = []
	procs = []
	with open(filename, "r") as f:
		for line in f.readlines():
			el = line.split()

			if line.find("[ntasks]") > -1:
				procs += [int(el[-1])]

			if len(el) == 6:
				if el[1] == "wallclock":
					totTime += [float(el[2])]
					wallTime += [float(el[3])]
				elif el[1] == "mpi":
					mpiTotTime += [float(el[2])]
					mpiTime += [float(el[3])]

	return procs, wallTime, totTime, mpiTime, mpiTotTime


def PlotTimings(filename):
	P, W, WT, M, MT = ParseIPMOutput(filename)

	fig = pylab.figure()
	ax = fig.add_subplot(111)
	ax2 = ax.twinx()

	ax.loglog(P, W, "b-", label="Wall clock")
	ax.loglog(P, W[1]*P[1]/array(P), "b--", label="Ideal wall clock")
	
	ax.loglog(P, M, "k-", label="MPI clock")
	ax2.semilogx(P, array(M)/array(W)*100, "k--", label="MPI clock %")

	ax.set_xlabel("Number of procs")
	ax.set_ylabel("Time [s]")
	ax2.set_ylabel("Time %")

	ax.set_xlim(1, P[0])

	#add merged legend
	lines = ax.get_lines()
	lines.extend(ax2.get_lines())
	labels = [l.get_label() for l in lines] 
	leg = ax.legend(lines, labels, loc="center left")

