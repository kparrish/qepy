from qepy import *

with pwx('./ex_si_q2r',
		title='si_pwx',
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
		) as pwCalc:
	with phx('./ex_si_q2r',
			title='si_ph',
			PWX=pwCalc,
			tr2_ph=1e-10,
			ldisp='.true.',
			nq1=2,
			nq2=2,
			nq3=2,
			amass=[1, 28],
			prefix='silicon',
			outdir='./tmp',
			fildyn='./si.dyn',
			walltime='24' # Can pass run parameters in the object creation
			) as phCalc:
		with q2r('./ex_si_q2r',
				title='si_q2r',
				PHX=phCalc,
				fildyn='./si.dyn',
				zasr='simple',
				flfrc='./harmonic_flfrc.dat'
				) as q2rCalc:
			try:
				q2rCalc.calculate(mode='queue', walltime='1')
			except (QepyRunning, QepySubmitted):
				pass


