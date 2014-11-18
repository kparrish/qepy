## ## ## qepy.py v.0.1.1
## ## ## Created: 06/12/2014 - KDP
import os
import commands
import numpy as np
from subprocess import Popen, PIPE
from shutil import rmtree
from os.path import isdir

import pwxlist as pwxl			# parameter lists of pwx
from qepyrc import *			# configuration file
from qepy_exceptions import *	# exceptions definitions


class pwx:
	"""
	qepy.pwx()
	A Quantum Espresso pw.x instance
	"""
	def __init__(self, qedir, **kwargs):
		usrHome = os.path.expanduser('~')
		self.qedir = qedir.replace('~', usrHome)
		self.cwd = os.getcwd()
		self.quote_control_params = {}
		self.control_params = {}
		self.paren_system_params = {}
		self.system_params = {}
		self.electrons_params = {}
		self.paren_electrons_params = {}
		self.ions_params = {}
		self.cell_params = {}
		self.atomic_species_params = {}
		self.atomic_positions_params = {}
		self.k_points_params = {}
		self.cell_parameters_params = {}
		self.constraints_params = {}
		self.occupations_params = {}
		self.atomic_forces_params = {}
		for key in pwxl.quote_control_keys:
			self.quote_control_params[key] = None
		for key in pwxl.control_keys:
			self.control_params[key] = None
		for key in pwxl.paren_system_keys:
			self.paren_system_params[key] = None
		for key in pwxl.system_keys:
			self.system_params[key] = None
		for key in pwxl.electrons_keys:
			self.electrons_params[key] = None
		for key in pwxl.paren_electrons_keys:
			self.paren_electrons_params[key] = None
		for key in pwxl.ions_keys:
			self.ions_params[key] = None
		for key in pwxl.cell_keys:
			self.cell_params[key] = None
		for key in pwxl.atomic_species_keys:
			self.atomic_species_params[key] = None
		for key in pwxl.atomic_positions_keys:
			self.atomic_positions_params[key] = None
		for key in pwxl.k_points_keys:
			self.k_points_params[key] = None
		for key in pwxl.cell_parameters_keys:
			self.cell_parameters_params[key] = None
		for key in pwxl.constraints_keys:
			self.constraints_params[key] = None
		for key in pwxl.occupations_keys:
			self.occupations_params[key] = None
		for key in pwxl.atomic_forces_keys:
			self.atomic_forces_params[key] = None

		self._set(**kwargs)

	def __enter__(self):
		"""
		On enter, if directory does not exist, create it.
		Change into directory.

		Syntax:
		with qepy() as calc:
			try:
				calc.my_command()
			except (QepyException):
				do something
		"""
		## Create Paths
		qedir = str(self.qedir).rstrip('/')
		outdir = self.quote_control_params['outdir'].strip(' \'\"\t\n\r/.')
		if not isdir(qedir):
			os.mkdir('{0}'.format(qedir))
			os.mkdir('{0}/{1}'.format(qedir, outdir)) 
		os.chdir(qedir)
		return self
	
	def __exit__(self, exc_type, exc_val, exc_tb):
		"""
		On exit, change back to original directory.
		"""
		os.chdir(self.cwd)
		return False	# Allows exception to propogate out

	def _set(self, **kwargs):
		"""
		Set kwargs parameters to those defined by pwxlist.
		"""
		for key in kwargs:
			if self.quote_control_params.has_key(key):
				self.quote_control_params[key] = '\'{0}\''.format(kwargs[key].strip('\'\" '))
			elif self.control_params.has_key(key):
				self.control_params[key] = kwargs[key]
			elif self.paren_system_params.has_key(key):
				self.paren_system_params[key] = kwargs[key]
			elif self.system_params.has_key(key):
				self.system_params[key] = kwargs[key]
			elif self.electrons_params.has_key(key):
				self.electrons_params[key] = kwargs[key]
			elif self.paren_electrons_params.has_key(key):
				self.paren_electrons_params[key] = kwargs[key]
			elif self.ions_params.has_key(key):
				self.ions_params[key] = kwargs[key]
			elif self.cell_params.has_key(key):
				self.cell_params[key] = kwargs[key]
			elif self.atomic_species_params.has_key(key):
				self.atomic_species_params[key] = kwargs[key]
			elif self.atomic_positions_params.has_key(key):
				self.atomic_positions_params[key] = kwargs[key]
			elif self.k_points_params.has_key(key):
				self.k_points_params[key] = kwargs[key]
			elif self.cell_parameters_params.has_key(key):
				self.cell_parameters_params[key] = kwargs[key]
			elif self.constraints_params.has_key(key):
				self.constraints_params[key] = kwargs[key]
			elif self.occupations_params.has_key(key):
				self.occupations_params[key] = kwargs[key]
			elif self.atomic_forces_params.has_key(key):
				self.atomic_forces_params[key] = kwargs[key]
			else:
				raise TypeError('Parameter not defined: '+ key)

	def calculate(self, **kwargs):
		inFile = self.quote_control_params['title'].strip('\'\"') + '.in'
		outFile = self.quote_control_params['title'].strip('\'\"') + '.out'
		if not os.path.isfile(inFile):
			## Create input file
			_qeControl(self)
			_qeSystem(self)
			_qeElectrons(self)
			_qeIons(self)
			_qeAtomicSpecies(self)
			_qeAtomicPositions(self)
			_qeKpoints(self)
			_qeCellParameters(self)
			_qeOccupations(self)
			_qeAtomicForces(self)

		## Set run/parallelization commands
		# QEPYRC defaults
		self.run_params = {}
		self.parallel_params = {}
		for key in pwxl.run_keys:
			self.run_params[key] = None
		for key in pwxl.parallel_keys:
			self.parallel_params[key] = None
		for key in QEPYRC:
			if self.run_params.has_key(key):
				self.run_params[key] = QEPYRC[key]

		for key in kwargs:
			if self.run_params.has_key(key):
				self.run_params[key] = kwargs[key]
			elif self.parallel_params.has_key(key):
				self.parallel_params[key] = kwargs[key]
			else:
				raise TypeError('Parameter not defined: '+ key)
		if self.run_params['jobname'] == None:
			self.run_params['jobname'] = self.quote_control_params['title'].strip('\'\"')

		## Submit/Run job
		if self._job_in_queue():	# If in queue, exit
			raise QepyRunning()
		elif os.path.isfile(outFile):
			if self._energy() != False:	# If already done, exit
				pass
			else:
				raise QepyNotComplete('Job not found in queue, energy not found in out file')
		elif os.path.isfile('jobid'):	# jobid file already exists
			raise QepyNotComplete('Found jobid file, but not running and no output file')
		else:	# If not in queue and not done, run	
			if self.run_params['mode'] == 'queue':
				_qeSubmission(self)
				p = Popen(['qsub', 'pwxrun.sh'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
				out, err = p.communicate()
				if out == '' or err != '':
					raise Exception(err)
				else:
					f = open('jobid','w')
					f.write(out)
					f.close()
					raise QepySubmitted(out)
			elif self.run_params['mode'] == 'local':
				self._local(**kwargs)
	
	def _job_in_queue(self):
		if not os.path.isfile('jobid'):
			return False
		else:
			# Get jobid
			jobid = open('jobid').readline().strip()

			# See if jobid is in the queue
			jobids_in_queue = commands.getoutput('qselect').split('\n')
			if jobid in jobids_in_queue:
				status, output = commands.getstatusoutput('qstat {0}'.format(jobid))
				if status == 0:
					lines = output.split('\n')
					fields = lines[2].split()
					job_status = fields[4]
					if job_status == 'C':
						return False
					else:
						return True
			else:
				return False

	def _local(self, **kwargs):
		"""
		Run the calculation through command line
		"""
		inFileName = self.quote_control_params['title'].strip('\'\"') + '.in'
		outFileName = self.quote_control_params['title'].strip('\'\"') + '.out'
		command = self.run_params['pw.x']
		for key, value in self.parallel_params.items():
			if value is not None:
				command += ' -{0} {1} '.format(str(key), str(value))
		command += ' < {0} > {1}'.format(inFileName, outFileName)
		os.system(command)

	def _energy(self):
		fileName = self.quote_control_params['title'].strip('\'\"') + '.out'
		outFile = open(str(fileName), 'r')

		while True:
			myString = outFile.readline()

			if myString.find('!') > -1:
				outFile.close()
				return float(myString.split()[4])

			elif myString == '':
				outFile.close()
				break
		return False
	#-- END _energy --#

	def get_energy(self):
		# add if running, add if need to submit
		energy = self._energy()
		if not energy:
			raise QepyNotComplete()
		else:
			return energy

def _qeControl(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')

	inFile.write('# This file generated by qepy\n')
	inFile.write('# For more info, visit github.com/kparrish/qepy\n')
	inFile.write(' &control' + '\n')
	for key, val in self.quote_control_params.items():
		if val is not None:
			inFile.write('  {0:s}={1:s},\n'.format(str(key), str(val)))
	for key, val in self.control_params.items():
		if val is not None:
			inFile.write('  {0:s}={1:s},\n'.format(str(key), str(val)))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeControl --#

def _qeSystem(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')

	inFile.write(' &system' + '\n')
	for key, val in self.paren_system_params.items():
		if val is not None:
			if len(val) is not 2:
				raise ValueError('Value for {0} must be a list or tuple of length 2'.format(key))
			elif key is 'celldm':	# Have celldm default to an angstrom input
				inFile.write('  {0}({1})={2},\n'.format(str(key), str(val[0]), str(val[1]/0.5291772108)))
			else:
				inFile.write('  {0}({1})={2},\n'.format(str(key), str(val[0]), str(val[1])))
	for key, val in self.system_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeSystem --#

def _qeElectrons(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')

	inFile.write(' &electrons' + '\n')
	for key, val in self.electrons_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	for key, val in self.paren_electrons_params.items():
		if val is not None:
			if len(val) is not 2:
				raise ValueError('Value for {0} must be a list or tuple of length 2'.format(key))
			else:
				inFile.write('  {0}({1})={2},\n'.format(str(key), str(val[0]), str(val[1])))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeElectrons --#

def _qeIons(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')

	inFile.write(' &ions' + '\n')
	for key, val in self.ions_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeIons --#

def _qeAtomicSpecies(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')

	inFile.write('ATOMIC_SPECIES\n')
	asp = self.atomic_species_params['atomic_species']
	if not isinstance(asp[0], list):
		if len(asp) != 3:
			raise ValueError('\'atomic_species\' must be len 3 or list of items of len 3')
		else:
			inFile.write('  {0} {1} {2}\n'.format(str(asp[0]), str(asp[1]), str(asp[2])))
	else:
		for line in range(len(asp)):
			if len(asp[line]) != 3:
				raise ValueError('\'atomic_species\' must be len 3 or list of items of len 3')
			else:
				inFile.write('  {0} {1} {2}\n'.format(str(asp[line][0]), str(asp[line][1]), str(asp[line][2])))
	inFile.close()
#-- END _qeAtomicSpecies --#

def _qeAtomicPositions(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')
	
	ap = self.atomic_positions_params['atomic_positions']
	if ap.lower() != 'alat' and ap.lower() != 'bohr' and \
			ap.lower() != 'angstrom' and ap.lower() != 'crystal':
		raise ValueError('\'atomic_positions\' must be either \'alat\', \'bohr\', \'angstrom\', or \'crystal\'')
	inFile.write('ATOMIC_POSITIONS {0}\n'.format(str(ap)))
	apl = self.atomic_positions_params['atomic_positions_list']
	if not isinstance(apl[0], list):
		if len(apl) != 4 and len(apl) != 7:
			raise ValueError('\'atomic_species_list\' must be len 4 or 7 or list of items of len 4 or 7')
		elif len(apl) == 4:
			inFile.write('  {0} {1} {2} {3}\n'.format(apl[0], apl[1], apl[2], apl[3])) 
		else:
			inFile.write('  {0} {1} {2} {3} {4} {5} {6}\n'.format(apl[0], apl[1], apl[2], apl[3], apl[4], apl[5], apl[6])) 
	else:
		for line in range(len(apl)):
			if len(apl[line]) != 4 and len(apl[line]) !=7:
				raise ValueError('\'atomic_species_list\' must be len 4 or 7 or list of items of len 4 or 7')
			elif len(apl[line]) == 4:
				inFile.write('  {0} {1} {2} {3}\n'.format(apl[line][0], apl[line][1], apl[line][2], apl[line][3])) 
			else:
				inFile.write('  {0} {1} {2} {3} {4} {5} {6}\n'.format(apl[line][0], apl[line][1], apl[line][2], apl[line][3], apl[line][4], apl[line][5], apl[line][6])) 
	inFile.close()
#-- END _qeAtomicPositions --#

def _qeKpoints(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')

	kpt = self.k_points_params['k_points']
	inFile.write('K_POINTS {0}\n'.format(kpt))
	if kpt.lower() == 'tpiba' or kpt.lower() == 'crystal' or \
			kpt.lower() == 'tpiba_b' or kpt.lower() == 'crystal_b' or \
			kpt.lower() == 'tpiba_c' or kpt.lower() == 'crystal_c':
		kpl = self.k_points_params['k_points_list']
		inFile.write('  ')
		for line in range(len(kpl)):
			for item in range(len(kpl[line])):
				inFile.write('{0} '.format(kpl[line][item]))
			inFile.write('\n')
	elif kpt.lower() == 'automatic':
		kpl = self.k_points_params['k_points_list']
		inFile.write('  ')
		for item in range(len(kpl)):
			inFile.write('{0} '.format(kpl[item]))
		inFile.write('\n')
	elif kpt.lower() != 'gamma':
		raise ValueError('\'k_points\' must be either tpiba | crystal | tpiba_b | crystal_b | tipiba_c | crystal_c | automatic | gamma')
	inFile.close()
#-- END _qeKpoints --#

def _qeCellParameters(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')
	
	if self.cell_parameters_params['units'] != None:
		inFile.write('CELL_PARAMETERS {{ {0} }}\n'.format(self.cell_parameters_params['units']))
		if self.cell_parameters_params['v1'] is not None and \
				self.cell_parameters_params['v2'] is not None and \
				self.cell_parameters_params['v3'] is not None:
			v1 = self.cell_parameters_params['v1']
			v2 = self.cell_parameters_params['v2']
			v3 = self.cell_parameters_params['v3']
			if len(v1) != 3 or len(v2) != 3 or len(v3) != 3:
				raise ValueError('v1, v2, and v3 must have dimensions of 3')
			inFile.write(' {0:s} {1:s} {2:s}\n'.format(str(v1[0]), str(v1[1]), str(v1[2])))
			inFile.write(' {0:s} {1:s} {2:s}\n'.format(str(v2[0]), str(v2[1]), str(v2[2])))
			inFile.write(' {0:s} {1:s} {2:s}\n'.format(str(v3[0]), str(v3[1]), str(v3[2])))
	inFile.close()
#-- END qeCellParameters --#

def _qeConstraints(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')
	
	use = False
	for key, val in self.constraints_params.items():
		if val is not None:
			use = True
		
	if use:
		inFile.write('CONSTRAINTS\n')
		inFile.write('  {0}'.format(str(self.constraints_params['nconstr'])))
		if self.constraints_params['constr_tol'] is not None:
			inFile.write(' {0}'.format(str(self.constraints_params['constr_tol'])))
		inFile.write('\n')

		constr = self.constraints_params['constr']
		inFile.write('  ')
		if not isinstance(constr[0], list):
			for i in range(len(constr)):
				inFile.write('{0} '.format(str(constr[i])))
			inFile.write('\n')
		else:
			for i in range(len(constr)):
				for j in range(len(constr[i])):
					inFile.write('{0} '.format(str(constr[i][j])))
				inFile.write('\n')
	inFile.close()
#-- END _qeConstraints --#

def _qeOccupations(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')

	use = False
	for key, val in self.occupations_params.items():
		if val is not None:
			use = True
		
	if use:
		inFile.write('OCCUPATIONS\n')
		for key, val in self.occupations_params.items():
			if val is not None:
				if not isinstance(val[0], list):
					if len(val) > 10:
						raise ValueError(key + ' cannot must have len <= 10!')
					inFile.write('  ')
					for i in range(len(val)):
						inFile.write('{0} '.format(val[i]))
					inFile.write('\n')
				else:
					for i in range(len(val)):
						if len(val[i]) > 10:
							raise ValueError(key + ' nested list cannot must have len <= 10!')
						inFile.write('  ')
						for j in range(len(val[i])):
							inFile.write('{0} '.format(val[i][j]))
						inFile.write('\n')
	inFile.close()
#-- END _qeOccupations --#

def _qeAtomicForces(self):
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(str(fileName), 'a')
	
	use = False
	for key, val in self.atomic_forces_params.items():
		if val is not None:
			use = True
		
	if use:
		inFile.write('ATOMIC_FORCES\n')
		for key, val in self.atomic_forces_params.items():
			if val is not None:
				if not isinstance(val[0], list):
					if len(val) != 4:
						raise ValueError(key + ' must have len of 4!')
					inFile.write('  ')
					for i in range(len(val)):
						inFile.write('{0} '.format(val[i]))
					inFile.write('\n')
				else:
					for i in range(len(val)):
						if len(val[i]) != 4:
							raise ValueError(key + ' must have len of 4!')
						inFile.write('  ')
						for j in range(len(val[i])):
							inFile.write('{0} '.format(val[i][j]))
						inFile.write('\n')
	inFile.close()
#-- END _qeAtomicForces --#

def _qeSubmission(self):
	inFileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	outFileName = self.quote_control_params['title'].strip('\'\"') + '.out'
	runFile = open('pwxrun.sh', 'w')
	runFile.write('#!/bin/bash' + '\n')
	runFile.write('#PBS {0}\n'.format(str(self.run_params['options'])))
	runFile.write('#PBS -l nodes={0}:ppn={1}\n'.format(str(self.run_params['nodes']), self.run_params['ppn']))
	runFile.write('#PBS -l walltime={0}:00:00\n'.format(str(self.run_params['walltime'])))
	runFile.write('#PBS -l mem={0}\n'.format(str(self.run_params['mem'])))
	runFile.write('#PBS -l vmem={0}\n'.format(str(self.run_params['vmem'])))
	runFile.write('#PBS -N {0}\n'.format(str(self.run_params['jobname'])))
	runFile.write('cd $PBS_O_WORKDIR\n')
	runFile.write('\n')
	runFile.write('mpirun -np {0} {1} '.format(str(self.run_params['ppn']), str(self.run_params['pw.x'])))
	for key, value in self.parallel_params.items():
		if value is not None:
			runFile.write('-{0} {1} '.format(str(key), str(value)))
	runFile.write('< {0} > {1}\n'.format(str(inFileName), str(outFileName)))
	runFile.close()
#-- END _qeSubmission --#

