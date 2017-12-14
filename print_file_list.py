import os

for root, folders, files in os.walk('py_box3'):
	for file in files:
		if '.py' in file and '.pyc' not in file:
			print('{}\t{}'.format(os.path.join(root, file), file))