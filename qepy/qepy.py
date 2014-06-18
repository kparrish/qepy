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
		self.ions_params = {}
		self.cell_params = {}
		self.atomic_species_params = {}
		self.atomic_positions_params = {}
		self.k_points_params = {}
		self.cell_parameters_params = {}
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

		self.set(**kwargs)

	def set(self, **kwargs):
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
		_qeAtomicSpecies(self)
		_qeAtomicPositions(self)
		_qeKpoints(self)
		_qeCellParameters(self)
		# _qeSubmission(self)
		
	def energy(self):
		"""
		Reads the energy from Quantum Espresso output file.

		ntpy.dft.energy(fileName)
		Parameters
		----------
			fileName : str
				A string containing the QE filename to read
				energy from
		Returns
		-------
			energy : float
				The energy found in the QE file.
		"""
		dftFile = open(fileName)

		while True:
			myString = dfptFile.readline()

			if myString.find('!') > -1:
				dftFile.close()
				return float(myString.split()[4])

			# Add if job error found
			elif myString == '':
				dftFile.close()
				###################throw error here
				break
		return 0
	#-- END energy --#

def _qeControl(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"')
	inFile = open(qedir + '/' + str(fileName), 'a')
	inFile.write('# This file generated by qepy\n')
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
	fileName = self.quote_control_params['title'].strip('\'\"')
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
	fileName = self.quote_control_params['title'].strip('\'\"')
	inFile = open(qedir + '/' + str(fileName), 'a')
	inFile.write(' &electrons' + '\n')
	for key, val in self.electrons_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeElectrons --#

def _qeAtomicSpecies(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"')
	inFile = open(qedir + '/' + str(fileName), 'a')
	inFile.write('ATOMIC_SPECIES' + '\n')
	# inFile.write(' {0:s} {1:s} {2:s}\n'.format(atom_id, atom_mass, atom_pot))
	for key, val in self.atomic_species_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.close()
#-- END _qeAtomicSpecies --#

def _qeAtomicPositions(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"')
	inFile = open(qedir + '/' + str(fileName), 'a')
	inFile.write('ATOMIC_POSITIONS crystal' + '\n')
	# for i in range(numAtoms):
		# inFile.write(' {0:s} {1:s} {2:s} {3:s}\n'.format(\
			# str(atomSym[i]), str(atomPos[i,0]), str(atomPos[i,1]), str(atomPos[i,2])))
	for key, val in self.atomic_positions_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.close()
#-- END _qeAtomicPositions --#

def _qeKpoints(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.quote_control_params['title'].strip('\'\"')
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
	fileName = self.quote_control_params['title'].strip('\'\"')
	inFile = open(qedir + '/' + str(fileName), 'a')
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


