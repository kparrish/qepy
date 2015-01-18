from exceptions import Exception

class QepyException(Exception):
	def __init__(self, message=''):
		self.message = message
	
	def __str__(self):
		return self.message

class QepyNotComplete(QepyException):
	pass

class QepySubmitted(QepyException):
	def __init__(self, jobid):
		self.jobid = jobid
	
	def __str__(self):
		return repr(self.jobid)

class QepyRunning(QepyException):
	pass

class QepyCrash(QepyException):
	def __init__(self, message=''):
		self.message = message
	def __str__(self):
		return self.message

