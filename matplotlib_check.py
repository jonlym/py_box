import os

i = 0
for root, folders, files in os.walk('.'):
	for file in files:
		if '.py' in file and '.pyc' not in file:
			file_path = os.path.join(root, file)
			with open(file_path, 'r') as f_ptr:
				for line in f_ptr:
					if 'matplotlib' in line:
						print('{}\t{} imports matplotlib'.format(i, file_path))
						i += 1
						break