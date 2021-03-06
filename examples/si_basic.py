from qepy import *

with pwx('./ex_si_basic',
		title='si_basic',
		calculation='scf',
		prefix='silicon',
		outdir='./tmp',
		pseudo_dir='../',

		ibrav=2,
		celldm=[1, 5.431],
		nat=2,
		ntyp=1,
		ecutwfc=20.0,

		atomic_species=['Si', 28, 'Si.pz-vbc.UPF'],
		atomic_positions='alat',
		atomic_positions_list=[['Si', 0, 0, 0],
								['Si', 0.25, 0.25, 0.25]],
		k_points='automatic',
		k_points_list=[4, 4, 4, 1, 1, 1],
		) as calc:
	try:
		calc.calculate(mode='local')
	except (QepyRunning, QepySubmitted):
		pass
