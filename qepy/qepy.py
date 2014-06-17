## ## ## qepy.py v.0.0
## ## ## Created: 06/12/2014 - KDP
import os
import numpy as np
from shutil import rmtree
from os.path import isdir

control_keys = [
	'calculation',
	'title',
	'verbosity',
	'restart_mode',
	'wf_collect',
	'nstep',
	'iprint',
	'tstress',
	'tprnfor',
	'dt',
	'outdir',
	'wfcdir',
	'prefix',
	'lkpoint_dir',
	'max_seconds',
	'etot_conv_thr',
	'forc_conv_thr',
	'disk_io',
	'pseudo_dir',
	'tefield',
	'dipfield',
	'lelfield',
	'nberrycyc',
	'lorbm',
	'lberry',
	'gdir',
	'nppstr',
]

system_keys = [
	'ibrav',
	'celldm',	# parenthesis
	'nat',
	'ntyp',
	'nbnd',
	'tot_charge',
	'tot_magnetization',
	'ecutwfc',
	'ecutrho',
	'ecutfock',
	'nr1',
	'nr2',
	'nr3',
	'nr1s',
	'nr2s',
	'nr3s',
	'nosym',
	'nosym_evc',
	'noinv',
	'no_t_rev',
	'force_symmorphic',
	'use_all_frac',
	'occupations',
	'one_atom_occupations',
	'starting_spin_angle',
	'degauss',
	'smearing',
	'nspin',
	'noncolin',
	'ecfixed',
	'qcutz',
	'q2sigma',
	'input_dft',
	'exx_fraction',
	'screening_parameter',
	'exxdiv_treatment',
	'x_gamma_extrapolation',
	'ecutvcut',
	'nqx1',
	'nqx2',
	'nqx3',
	'lda_plus_u',
	'lda_plus_u_kind',
	'U_projection_type',
	'edir',
	'emaxpos',
	'eopreg',
	'eamp',
	'constrained_magnetization',
	'lambda',
	'report',
	'lspinorb',
	'assume_isolated',
	'esm_bc',
	'esm_w',
	'esm_efield',
	'esm_nfit',
	'vdw_corr',
	'london',
	'london_s6',
	'london_rcut',
	'xdm',
	'xdm_a1',
	'xdm_a2',
]

electrons_keys = [
	'electron_maxstep',
	'scf_must_converge',
	'conv_thr',
	'adaptive_thr',
	'conv_thr_init',
	'conv_thr_multi',
	'mixing_mode',
	'mixing_beta',
	'mixing_ndim',
	'mixing_fixed_ns',
	'diagonalization',
	'ortho_para',
	'diago_thr_init',
	'diago_cg_maxiter',
	'diago_david_ndim',
	'diago_full_acc',
	'efield',
	'startingpot',
	'startingwfc',
	'tqr',
]

ions_keys = [
	'ion_dynamics',
	'ion_positions',
	'phase_space',
	'pot_extrapolation',
	'wfc_extrapolation',
	'remove_rigid_rot',
	'ion_temperature',
	'tempw',
	'tolp',
	'delta_t',
	'nraise',
	'refold_pos'
	'upscale',
	'bfgs_ndim',
	'trust_radius_max',
	'trust_radius_min',
	'trust_radius_ini',
	'w_1',
	'w_2',
]

cell_keys = [
	'cell_dynamics',
	'press',
	'wmass',
	'cell_factor',
	'press_conv_thr',
	'cell_dofree',
]

atomic_species_keys = [
	'atom_symbols',
	'atom_masses',
	'pseudopot',
]

atomic_positions_keys = [
	# 'atom_symbols',
	'atom_positions',
	'if_pos',		#parenthesis
]


k_points_keys = [
	'K_POINTS',
	'nks',
	'xk_x',
	'xk_y',
	'xk_z',
	'wk',
	'nk1',
	'nk2',
	'nk3',
	'sk1',
	'sk2',
	'sk3',
]

cell_parameters_keys = [
	'v1',
	'v2',
	'v3',
]

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
		self.control_params = {}
		self.system_params = {}
		self.electrons_params = {}
		self.ions_params = {}
		self.cell_params = {}
		self.atomic_species_params = {}
		self.atomic_positions_params = {}
		self.k_points_params = {}
		self.cell_parameters_params = {}
		for key in control_keys:
			
			self.control_params[key] = None
		for key in system_keys:
			self.system_params[key] = None
		for key in electrons_keys:
			self.electrons_params[key] = None
		for key in ions_keys:
			self.ions_params[key] = None
		for key in cell_keys:
			self.cell_params[key] = None
		for key in atomic_species_keys:
			self.atomic_species_params[key] = None
		for key in atomic_positions_keys:
			self.atomic_positions_params[key] = None
		for key in k_points_keys:
			self.k_points_params[key] = None
		for key in cell_parameters_keys:
			self.cell_parameters_params[key] = None

		self.set(**kwargs)

	def set(self, **kwargs):
		for key in kwargs:
			if self.control_params.has_key(key):
				self.control_params[key] = kwargs[key]
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
		outdir = self.control_params['outdir'].strip(' \t\n\r/.')
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
	fileName = self.control_params['title']
	inFile = open(qedir + '/' + str(fileName), 'a')
	inFile.write('# This file generated by qepy\n')
	inFile.write(' &control' + '\n')
	for key, val in self.control_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeControl --#

def _qeSystem(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.control_params['title']
	inFile = open(qedir + '/' + str(fileName), 'a')
	inFile.write(' &system' + '\n')
	for key, val in self.system_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeSystem --#

def _qeElectrons(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.control_params['title']
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
	fileName = self.control_params['title']
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
	fileName = self.control_params['title']
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
	fileName = self.control_params['title']
	inFile = open(qedir + '/' + str(fileName), 'a')
	# inFile.write('K_POINTS automatic' + '\n')
	# inFile.write(' {0} {0} {0} 0 0 0' + '\n'.format(str(kpt)))
	for key, val in self.k_points_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
	inFile.close()
#-- END _qeKpoints --#

def _qeCellParameters(self):
	qedir = str(self.qedir).rstrip('/')
	fileName = self.control_params['title']
	inFile = open(qedir + '/' + str(fileName), 'a')
	inFile.write('CELL_PARAMETERS cubic' + '\n')
	# for i in range(latVec[:,0].size):
		# inFile.write(' {0:s} {1:s} {2:s}\n'.format(str(latVec[i,0]), str(latVec[i,1]), str(latVec[i,2])))
	for key, val in self.cell_parameters_params.items():
		if val is not None:
			inFile.write('  {0}={1},\n'.format(str(key), str(val)))
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


