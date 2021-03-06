import os

QEPYRC = {'pw.x':'/opt/espresso-5.1/bin/pw.x',
		'ph.x':'/opt/espresso-5.1/bin/ph.x',
		'q2r.x':'/opt/espresso-5.1/bin/q2r.x',
		'mode':'queue',
		'command':'qsub',
		'options':'-joe',
		'walltime':'24',
		'nodes':1,
		'ppn':1,
		'mem':'2GB',
		'vmem':'2GB',
		'jobname:':'None'}

def read_configuration(fname):
	f = open(fname)
	for line in f:
		line = line.strip()

		if line.startswith('#'):
			pass
		elif line == '':
			pass
		else:
			if '#' in line:
				# take the part before the first #
				line = line.split('#')[0]
			key, value = line.split('=')
			QEPYRC[key.strip()] = value.strip()

config_files = [os.path.join(os.environ['HOME'], '.qepyrc'),
				'.qepyrc']

for cf in config_files:
	if os.path.exists(cf):
		read_configuration(cf)
