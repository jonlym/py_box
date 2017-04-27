# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 11:24:22 2017

@author: Jonathan Lym
"""

from warnings import warn
from os import getcwd, chdir
from os.path import isdir, exists

def check_files(job_file = 'joblist.txt', folder_file = './folderlist.txt'):
    home_dir = getcwd()
    with open(job_file, 'r') as job_ptr, open(folder_file, 'r') as folder_ptr:
        for job, folder in zip(job_ptr, folder_ptr):
            #Check if folder exists
            job = job.replace('\n', '')
            folder = folder.replace('\n', '')
            folder_exists = isdir(folder)
            if folder_exists:
                chdir(folder)
                job_exists = exists('%s.py' % job)
                if job_exists:
                    print 'Job %s.py at %s exists!' % (job, folder)
                else:
                    print 'Warning: Job %s.py does not exist at %s.' % (job, folder)
            else:
                print 'Warning: Folder %s does not exist.' % folder
            chdir(home_dir)
