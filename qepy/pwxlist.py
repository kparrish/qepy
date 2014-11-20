control_keys = [
	'wf_collect',
	'iprint',
	'tstress',
	'tprnfor',
	'dt',
	'lkpoint_dir',
	'max_seconds',
	'etot_conv_thr',
	'forc_conv_thr',
	'tefield',
	'dipfield',
	'lelfield',
	'nberrycyc',
	'lorbm',
	'lberry',
	'gdir',
	'nppstr',
]

quote_control_keys = [
	'calculation',
	'title',
	'verbosity',
	'restart_mode',
	'nstep',
	'outdir',
	'wfcdir',
	'prefix',
	'disk_io',
	'pseudo_dir',
]

system_keys = [
	'ibrav',
	'A',
	'B',
	'C',
	'cosAB',
	'cosAC',
	'cosBC',
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

paren_system_keys = [
	'celldm',
	'starting_magnetization',
	'Hubbard_U',
	'Hubbard_J0',
	'Hubbard_alpha',
	'Hubbard_beta',
	'Hubbard_J', #need to confirm
	'starting_ns_eigenvalue', #need to confirm
	'angle1',
	'angle2',
	'fixed_magnetization',
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

paren_electrons_keys = [
	'efield_cart',
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
	'atomic_species',
]

atomic_positions_keys = [
	'atomic_positions',
	'atomic_positions_list',
]

k_points_keys = [
	'k_points',
	'k_points_list',
]

cell_parameters_keys = [
	'units',
	'v1',
	'v2',
	'v3',
]

occupations_keys = [
	'f_inp1',
	'f_inp2',
]

constraints_keys = [
	'nconst',
	'constr_tol',
	'constr',
]

atomic_forces_keys = [
	'atomic_forces',
]

