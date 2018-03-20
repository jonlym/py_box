import subprocess


def parse_qstat():
	lines = subprocess.check_output(['qstat']).split('\n')
	for line in lines:
		if '-----' in line or line[0].isalpha():
			continue
		fields = [field for field in line.split(' ') if line != '']
		