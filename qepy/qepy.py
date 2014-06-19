## ## ## qepy.py v.0.0
## ## ## Created: 06/12/2014 - KDP
import os
import numpy as np
from shutil import rmtree
from os.path import isdir
import pwscflist as pwscfl

class PWscf:
	"""
	qepy.PWscf()
	Parameters
	----------
		restart_mode : str
			Switch
	Returns
	-------
		energy : float
			The energy found in the QE file.
	"""
	def __init__(self, qedir, **kwargs):
		usrHome = os.path.expanduser('~')
		self.qedir = qedir.replace('~', usrHome)
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
		for key in pwscfl.quote_control_keys:
			self.quote_control_params[key] = None
		for key in pwscfl.control_keys:
			self.control_params[key] = None
		for key in pwscfl.paren_system_keys:
			self.paren_system_params[key] = None
		for key in pwscfl.system_keys:
			self.system_params[key] = None
		for key in pwscfl.electrons_keys:
			self.electrons_params[key] = None
		for key in pwscfl.paren_electrons_keys:
			self.paren_electrons_params[key] = None
		for key in pwscfl.ions_keys:
			self.ions_params[key] = None
		for key in pwscfl.cell_keys:
			self.cell_params[key] = None
		for key in pwscfl.atomic_species_keys:
			self.atomic_species_params[key] = None
		for key in pwscfl.atomic_positions_keys:
			self.atomic_positions_params[key] = None
		for key in pwscfl.k_points_keys:
			self.k_points_params[key] = None
		for key in pwscfl.cell_parameters_keys:
			self.cell_parameters_params[key] = None
		for key in pwscfl.constraints_parameters_keys:
			self.constraints_params[key] = None
		for key in pwscfl.occupations_parameters_keys:
			self.occupations_params[key] = None
		for key in pwscfl.atomic_forces_parameters_keys:
			self.atomic_forces_params[key] = None

		self._set(**kwargs)

	def _set(self, **kwargs):
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

	def calculate(self, recalc=False):
		qedir = str(self.qedir).rstrip('/')
		outdir = self.quote_control_params['outdir'].strip(' \'\"\t\n\r/.')
		if not recalc:
			if isdir(qedir):
				raise OSError('Directory {0} exists, set recalc=False to override and delete the directory'.format(qedir))
		else:
			if isdir(qedir):
				rmtree(qedir)

		os.mkdir('{0}'.format(qedir))
		os.mkdir('{0}/{1}'.format(qedir, outdir)) #need to handle relative and absolute

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
		# _qeSubmission(self)
		
	def _energy(self):
		fileName = self.quote_control_params['title'].strip('\'\"') + '.out'
		outFile = open(qedir + '/' + str(fileName), 'r')

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
		energy = _energy(self)
		if not energy:
			raise QepyNotComplete('')
		else:
			return energy

def _qeControl(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

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
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

	inFile.write(' &system' + '\n')
	for key, val in self.paren_system_params.items():
		if val is not None:
			if len(val) is not 2:
				raise ValueError('Value for {0} must be a list or tuple of length 2'.format(key))
			else:
				inFile.write('  {0}({1})={2},\n'.format(str(key), str(val[0]), str(val[1])))
	for key, val in self.system_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeSystem --#

def _qeElectrons(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

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
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

	inFile.write(' &ions' + '\n')
	for key, val in self.ions_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeIons --#

def _qeAtomicSpecies(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

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
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')
	
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
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

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
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')
	
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
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

	
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
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

	
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
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"') + '.in'
	inFile = open(qedir + '/' + str(fileName), 'a')

	
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
	runFile = open(str(fileName), 'w')
	runFile.write('#!/bin/bash' + '\n')
	runFile.write('#PBS -j oe' + '\n')
	runFile.write('#PBS -l nodes=1:ppn={0:s}\n'.format(str(numProcs)))
	runFile.write('#PBS -l walltime={0:s}:00:00' + '\n'.format(str(wallHrs)))
	runFile.write('#PBS -l mem={0:s}gb' + '\n'.format(str(mem)))
	runFile.write('#PBS -l vmem={0:s}gb' + '\n'.format(str(vmem)))
	runFile.write('#PBS -N {0:s}\n'.format(str(title)))
	runFile.write('cd $PBS_O_WORKDIR' + '\n')
	runFile.write('\n')
	runFile.close()
#-- END _qeSubmission --#


