"""
YieldCondition
==============

"""

class YieldConditionBase(object):
	def __init__(self):
		raise NotImplementedError("Implement on deriving class")
	
	def __call__(self, prop):
		raise NotImplementedError("Implement on deriving class")
	
	
class YieldConditionConstantNormDecay(YieldConditionBase):
	"""	Yield True if abs(||psi(t)|| - ||psi(t+dt)||) > tol for given tolerance 
		
	"""
	
	def __init__(self, initialNorm, normDecayTolerance):
		self.OldNorm = initialNorm
		self.NormDecayTolerance = normDecayTolerance
		
	def __call__(self, prop):
		curNorm = prop.psi.GetNorm()
		normDiff = abs(curNorm - self.OldNorm)
		self.OldNorm = curNorm
		
		return normDiff > self.NormDecayTolerance