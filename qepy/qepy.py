## ## ## qepy.py v.0.0
## ## ## Created: 06/12/2014 - KDP
import os
import numpy as np

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
	'celldm(1)',
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
	'if_pos(1)',
	'if_pos(2)',
	'if_pos(3)',
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
		self.qedir = qedir
		self.calculator_params = {}
		self.system_params = {}
		self.electrons_params = {}
		self.ions_params = {}
		self.cell_params = {}
		self.atomic_species_params = {}
		self.atomic_positions_params = {}




	
	# Cleanup inputs
	if outdir == None:
		outdir=str(title)+'/tmp/'
	if prefix == None:
		prefix=str(title)

	numProcs = str(numProcs)
	wallHrs = str(wallHrs)
	mem = str(mem)
	vmem = str(vmem)
	npool = str(npool)
	kpt = str(kpt)
	ecutwfc = str(ecutwfc)
	atom_id = str(atom_id)
	atom_mass = str(atom_mass)
	atom_pot = str(atom_pot)

	# Read structure directory
	dfptFile = open(struct_dir + '/dfpt')

	latVec = np.zeros((3,3))
	latVecNum = -1
	latConst = 0.0
	numAtoms = -1
	atomSym = []
	qPoints = np.zeros((3))

	i = 0
	while True:
		myString = dfptFile.readline()

		if i == 0:
			latConst = float(myString)
			i += 1

		if myString.find('Lattice Vectors:') > -1:
			lineCont = dfptFile.readline().split()
			latVec[0,:] = float(lineCont[0]), float(lineCont[1]), float(lineCont[2])
			i += 1
			lineCont = dfptFile.readline().split()
			latVec[1,:] = float(lineCont[0]), float(lineCont[1]), float(lineCont[2])
			i += 1
			lineCont = dfptFile.readline().split()
			latVec[2,:] = float(lineCont[0]), float(lineCont[1]), float(lineCont[2])
			i += 1

		if myString.find('Number of Atoms:') > -1:
			lineCont = myString.split()
			numAtoms = int(lineCont[3])
			atomPos = np.zeros((numAtoms, 3))
			for atom in range(numAtoms):
				lineCont = dfptFile.readline().split()
				atomSym.append(lineCont[0])
				atomPos[atom, :] = lineCont[1:4]
				i += 1

		if myString.find('q Points needed') > -1:
			lineCont = dfptFile.readline().split()
			qPoints[:] = lineCont[:]
			i += 1

		if myString == '':
			break
	
	dfptFile.close()

	# Make executable directories
	if os.path.isdir(title):
		os.system('rm -r '+ title)
	if os.path.isdir(outdir):
		os.system('rm -r '+ outdir)
	os.system('mkdir ' + title)
	os.system('mkdir ' + outdir)
	
	# Make input file
	celldm = str(float(latConst) * float(ascale) / (0.5291772108))
	inFileName = './' + title+'/qe.in'

	_qeControl(inFileName,
			title=title, 
			calculation=calculation,
			restart_mode=restart_mode,
			outdir=outdir,
			prefix=prefix,
			tprnfor=tprnfor,
			etot_conv_thr=etot_conv_thr)

	_qeSystem(inFileName,
			ibrav=ibrav, 
			celldm=celldm,
			numAtoms=numAtoms,
			ntyp=ntyp,
			ecutwfc=ecutwfc)

	_qeElectrons(inFileName,
			conv_thr=conv_thr,
			mixing_beta=mixing_beta)

	_qeAtomicSpecies(inFileName,
			atom_id=atom_id,
			atom_mass=atom_mass,
			atom_pot=atom_pot)

	_qeAtomicPositions(inFileName,
			numAtoms=numAtoms,
			atomSym=atomSym,
			atomPos=atomPos)

	_qeKpoints(inFileName,
			kpt=kpt)

	_qeCellParameters(inFileName,
			latVec=latVec)

	# Make submission file
	runFileName = './' + title + '/run.sh'

	_qeSubmission(runFileName,
			title=title,
			numProcs=numProcs,
			wallHrs=wallHrs,
			mem=mem,
			vmem=vmem)
#-- END converge --#
		
def _qeControl(fileName,
			title=None, 
			calculation=None,
			restart_mode=None,
			outdir=None,
			prefix=None,
			tprnfor=None,
			etot_conv_thr=None):

	inFile = open(str(fileName), 'a')
	inFile.write(' &control' + '\n')
	inFile.write('  ' + 'title=\'{0:s}\'\n'.format(title))
	inFile.write('  ' + 'calculation=\'{0:s}\'\n'.format(calculation))
	inFile.write('  ' + 'restart_mode=\'{0:s}\'\n'.format(restart_mode))
	inFile.write('  ' + 'outdir=\'{0:s}\'\n'.format(outdir))
	# inFile.write('  ' + 'pseudo_dir=\'{0:s}\'\n'.format(pseudo_dir))
	inFile.write('  ' + 'prefix=\'{0:s}\'\n'.format(prefix))
	inFile.write('  ' + 'tprnfor={0:s}\n'.format(tprnfor))
	inFile.write('  ' + 'etot_conv_thr={0:s}\n'.format(etot_conv_thr))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeControl --#

def _qeSystem(fileName,
			ibrav=None, 
			celldm=None,
			numAtoms=None,
			ntyp=None,
			ecutwfc=None):

	inFile = open(str(fileName), 'a')
	inFile.write(' &system' + '\n')
	inFile.write('  ' + 'ibrav={0:s}\n'.format(ibrav))
	inFile.write('  ' + 'celldm(1)={0:s}\n'.format(celldm))
	inFile.write('  ' + 'nat={0:s}\n'.format(str(numAtoms)))
	inFile.write('  ' + 'ntyp={0:s}\n'.format(ntyp))
	inFile.write('  ' + 'ecutwfc={0:s}\n'.format(ecutwfc))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeSystem --#

def _qeElectrons(fileName,
			conv_thr=None,
			mixing_beta=None):

	inFile = open(str(fileName), 'a')
	inFile.write(' &electrons' + '\n')
	inFile.write('  ' + 'conv_thr={0:s}\n'.format(conv_thr))
	inFile.write('  ' + 'mixing_beta={0:s}\n'.format(mixing_beta))
	inFile.write(' /' + '\n')
	inFile.close()
#-- END _qeElectrons --#

def _qeAtomicSpecies(fileName,
			atom_id=None,
			atom_mass=None,
			atom_pot=None):

	inFile = open(str(fileName), 'a')
	inFile.write('ATOMIC_SPECIES' + '\n')
	inFile.write(' {0:s} {1:s} {2:s}\n'.format(atom_id, atom_mass, atom_pot))
	inFile.close()
#-- END _qeAtomicSpecies --#

def _qeAtomicPositions(fileName,
			numAtoms=None,
			atomSym=None,
			atomPos=None):

	inFile = open(str(fileName), 'a')
	inFile.write('ATOMIC_POSITIONS crystal' + '\n')
	for i in range(numAtoms):
		inFile.write(' {0:s} {1:s} {2:s} {3:s}\n'.format(\
			str(atomSym[i]), str(atomPos[i,0]), str(atomPos[i,1]), str(atomPos[i,2])))
	inFile.close()
#-- END _qeAtomicPositions --#

def _qeKpoints(fileName,
			kpt=None):

	inFile = open(str(fileName), 'a')
	inFile.write('K_POINTS automatic' + '\n')
	inFile.write(' {0} {0} {0} 0 0 0' + '\n'.format(str(kpt)))
	inFile.close()
#-- END _qeKpoints --#

def _qeCellParameters(fileName,
			latVec=None):

	inFile = open(str(fileName), 'a')
	inFile.write('CELL_PARAMETERS cubic' + '\n')
	for i in range(latVec[:,0].size):
		inFile.write(' {0:s} {1:s} {2:s}\n'.format(str(latVec[i,0]), str(latVec[i,1]), str(latVec[i,2])))
	inFile.close()
#-- END qeCellParameters --#

def _qeSubmission(fileName,
			title=None,
			numProcs=None,
			wallHrs=None,
			mem=None,
			vmem=None):

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

def _qeNumDisp():
	ndisp = 0
	while(True):
		fname = './structure/disp_{0}.xyz'.format(str(ndisp))
		fileExists = os.path.isfile(fname)
		if fileExists:
			ndisp += 1
		else:
			break
	return ndisp
#-- END _qeNumDisp --#

def _qeAcell():
	fname = './structure/acell'
	latVec = np.zeros((3,3))
	acell = open(fname)
	for i in range(3):
		lineCont = acell.readline().split
		latVec[0,:] = float(lineCont[0]), float(lineCont[1]), float(lineCont[2])
	return latVec
#-- END _qeAcell --#

def _qeNumAtomsDisp():
	natom = 0
	fname = './structure/disp_0.xyz'
	disp = open(fname)
	lineCont = acell.readline().split
	natom = int(lineCont[0])
	disp.close()
	return natom
#-- END _qeNumAtomsDisp --#

def _qeDisp(ndisp, natom):
	pos = np.zeros((3, natom, ndisp))
	symbol = np.empty((natom, ndisp), dtype=object)
	for d in range(ndisp): # displacement file
		fname = './structure/disp_{0}.xyz'.format(str(d))
		disp = open(fname)
		lineCont = disp.readline() # read two header lines
		lineCont = disp.readline()
		for a in range(natom):
			lineCont = disp.readline().split
			symbol[a,d] = lineCont[0]
			pos[:,a,d] = float(lineCont[2]), float(lineCont[3]), float(lineCont[4]) 
		disp.close()
	return pos, symbol
#-- END _qeDisp --#


def energy(fileName):
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

		if myString == '':
			print 'Did not find total energy!'
			dftFile.close()
			break
	return 0
#-- END energy --#
