from qepy import *

kpts = [2, 4, 6]
for i, k in enumerate(kpts):
	with pwx('./ex_si_kpts_' + str(k),
					title='si_kpts_' + str(k),
					calculation='scf',
					prefix='si_kpts_' + str(k),
					outdir='./tmp',

					ibrav=2,
					celldm=[1, 10.2],
					nat=2,
					ntyp=1,
					ecutwfc=20.0,

					atomic_species=['Si', 28, 'Si.pz-vbc.UPF'],
					atomic_positions='alat',
					atomic_positions_list=[['Si', 0, 0, 0],
										['Si', 0.25, 0.25, 0.25]],
					k_points='automatic',
					k_points_list=[k, k, k, 1, 1, 1],
				) as calc:
		try:
			calc.calculate(recalc=True)
		except(QepyException):
			pass

