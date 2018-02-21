#!/usr/bin/env python

import os
import argparse
from ase.io import read

parser = argparse.ArgumentParser(description = "Iterates through folderlist to find OUTCAR energies.")
parser.add_argument("folder_file", nargs = '?' default = './folderlist.txt', help = "File that contains the folder list of job array. Default is ./folderlist.txt")
args = parser.parse_args()

with open(args.folderlist, 'r') as f_ptr:
    for line in f_ptr:
        line = line.replace('\n', '')
        try:
            atoms = read(os.path.join(*['.', line, 'OUTCAR']))
        except:
            print 'Unable to process {}/OUTCAR'.format(line)
        else:
            print '{}\t{}'.format(line, atoms.get_potential_energy())
