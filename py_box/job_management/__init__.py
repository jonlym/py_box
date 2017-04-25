# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 21:26:06 2017

@author: Jon Lym
"""

from subprocess import call

parameters_dict = {'Farber': {'n_core': 20},
                   'Squidward': {'n_core': 16}}

class job_management(object):
    def __init__(self, jobs = []):
        self._jobs = []
        for job in jobs:
            self.append(job)
            
    def append(self, new_job):
        self._jobs.append(new_job)

    def extend(self, new_jobs):
        self._jobs.extend(new_jobs)

    def index(self, id_key = None, name_key = None):
        for i, job in enumerate(self._jobs):
            if id_key is not None:
                if job.job_id == id_key:
                    return i
            elif name_key is not None:
                if job.name == name_key:
                    return i

    def remove(self, id_key = None, name_key = None):
        for i, job in enumerate(self._jobs):
            if id_key is not None:
                if job.job_id == id_key:
                    self._jobs.pop(i)
            elif name_key is not None:
                if job.name == name_key:
                    self._jobs.pop(i)

    def __len__(self):
        return len(self._jobs)
    
    def __setitem__(self, index, job):
        self._jobs[index] = job

    def __getitem__(self, index):
        return self._jobs[index]
        
#    def __str__(self):
#        symbols = ''
#        for thermdat in self:
#            symbols += '%s, ' % thermdat.symbol
#        symbols = symbols[:-2]
#        return symbols

    
class job(object):
    def __init__(self, job_id = None, path = None, file_name = None, status = None, history = None):
        self.set_job_id(job_id)
        self.set_path(path)
        self.set_file_name(file_name)
        self.set_status(status)
        self.set_history(history)
        
    def set_job_id(self, job_id):
        self._job_id = job_id
        
    def set_path(self, path):
        self._path = path
        
    def set_file_name(self, file_name):
        self._file_name= file_name
        
    def set_status(self, status):
        self._status = status
        
    def set_history(self, history):
        self._history = history

    def get_job_id(self):
        return self._job_id
        
    def get_path(self):
        return self._path
        
    def get_file_name(self, file_name):
        return self._file_name
        
    def get_status(self):
        return self._status
        
    def get_history(self):
        return self._history
 
        
def read_qstat(job_management):
    """Calls the command qstat and reads the job listings."""
    out = subprocess(['qstat', '-u', 'jlym'])
    lines = out.split('\n')
    for line in lines:
        job_attr = line.split(' ')
        job_attr = [x for x in job_attr if x != '']
        Job = job(job_id = job_attr[0], status = job_attr[4])
        job_management.append(Job)
        
def read_log_file(log_path):
    with open(log_path, 'r') as log_file:
        
        