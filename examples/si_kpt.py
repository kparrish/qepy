from qepy import *
import numpy as np
import matplotlib.pyplot as plt

kpts = [2, 4, 6, 8, 10]
energy = np.zeros((len(kpts)))

for i, k in enumerate(kpts):
	with pwx('./ex_si_kpts_' + str(k),
					title='si_kpts_' + str(k),
					calculation='scf',
					prefix='si_kpts_' + str(k),
					outdir='./tmp',
					pseudo_dir='../',

					ibrav=2,
					celldm=[1, 5.431],
					nat=2,
					ntyp=1,
					ecutwfc=20.0,

					atomic_species=['Si', 28, 'Si.pz-bhs.UPF'],
					atomic_positions='alat',
					atomic_positions_list=[['Si', 0, 0, 0],
										['Si', 0.25, 0.25, 0.25]],
					k_points='automatic',
					k_points_list=[k, k, k, 1, 1, 1],
				) as calc:
		try:
			calc.calculate()
			energy[i] = calc.get_energy()
		except(QepyException):
			pass

plt.scatter(kpts, energy)
plt.xlabel('Number of kpts')
plt.ylabel('Total energy')
plt.savefig('pic_kpt.png')
plt.clf()
