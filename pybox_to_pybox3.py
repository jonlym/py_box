import os

for root, folders, files in os.walk('./py_box3'):
	for file in files:
		if '.py' in file and '.pyc' not in file:
			print('Processing {}'.format(file))
			file_path = os.path.join(root, file)
			lines = []
			with open(file_path, 'r') as f_ptr:
				for line in f_ptr:
					lines.append(line.replace('from py_box ', 'from py_box3 '))
			with open(file_path, 'w') as f_ptr:
				f_ptr.write(''.join(lines))