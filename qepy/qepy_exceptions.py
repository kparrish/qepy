from exceptions import Exception

class QepyException(Exception):
	pass

class QepyNotComplete(Exception):
	def __init__(self, message=''):
		self.message = message
	
	def __str__(self):
		return self.message
