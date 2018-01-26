# -*- coding: utf-8 -*-
"""
Created on Wed May 03 11:35:25 2017

@author: Jonathan Lym
"""

import os, sys
import numpy as np

from ase.utils import devnull, str
from ase.calculators.vasp import Vasp

class Vasp(Vasp):
    def run(self):
        if Vasp.track_output:
            Vasp.out = Vasp.output_template + str(self.run_counts) + '.out'
            Vasp.run_counts += 1
        else:
            Vasp.out = Vasp.output_template + '.out'
        stderr = sys.stderr
        p = self.input_params
        if p['txt'] is None:
            sys.stderr = devnull
        elif p['txt'] == '-':
            pass
        elif isinstance(p['txt'], str):
            sys.stderr = open(p['txt'], 'w')
        if 'VASP_COMMAND' in os.environ:
            vasp = os.environ['VASP_COMMAND']
            exitcode = os.system('%s >> %s' % (vasp, Vasp.out))
        elif 'VASP_SCRIPT' in os.environ:
            vasp = os.environ['VASP_SCRIPT']
            locals = {}
            exec(compile(open(vasp).read(), vasp, 'exec'), {}, locals)
            exitcode = locals['exitcode']
        else:
            raise RuntimeError('Please set either VASP_COMMAND'
                               ' or VASP_SCRIPT environment variable')
        sys.stderr = stderr
        if exitcode != 0:
            raise RuntimeError('Vasp exited with exit code: %d.  ' % exitcode)