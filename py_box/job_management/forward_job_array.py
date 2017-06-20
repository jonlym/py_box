#!/usr/bin/env python
"""
Created on Thu Apr 27 11:24:22 2017

@author: Jonathan Lym
"""

import argparse
from warnings import warn
from os import getcwd, chdir
from os.path import isdir, exists

#Parse arguments
parser = argparse.ArgumentParser(description = "Checks the paths in specified by job_in and folder_in and writes the unfinished jobs to job_out and folder_out.")
parser.add_argument("job_in", default = './joblist.txt', help = "File that holds the input jobs.")
parser.add_argument("folder_in", default = './folderlist.txt', help = "File that contains paths to input jobs.")
parser.add_argument("job_out", default = './joblist_out.txt', help = "File that uncompleted jobs will be written to.")
parser.add_argument("folder_out", default = './folderlist_out.txt', help = "File that paths to uncomplete jobs will be written to.")
args = parser.parse_args()
job_in = args.job_in
folder_in = args.folder_in
job_out = args.job_out
folder_out = args.folder_out

i = 0
home_dir = getcwd()
with open(job_in, 'r') as job_in_ptr, \
        open(folder_in, 'r') as folder_in_ptr,\
        open(job_out, 'w') as job_out_ptr,\
        open(folder_out, 'w') as folder_out_ptr:
    for job, folder in zip(job_in_ptr, folder_in_ptr):
        #Check if folder exists
        job = job.replace('\n', '')
        folder = folder.replace('\n', '')
        folder_exists = isdir(folder)
        if folder_exists:
            chdir(folder)
            job_exists = exists('{}.out'.format(job))
            if job_exists:
                if 'Completed {}.py'.format(job) not in open('{}.out'.format(job), 'r').read():
                    print '{}\t{}.py did not complete. Adding to job_out paths.'.format(i, job)
                    job_out_ptr.write('{}\n'.format(job))
                    folder_out_ptr.write('{}\n'.format(folder))
                    i += 1
            else:
                print '{}.out does not exist at {}'.format(job, folder)
        chdir(home_dir)
