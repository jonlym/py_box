# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:38:29 2016

@author: Jonathan Lym
"""

from ase.calculators.vasp import Vasp

def set_geo_calc_ZnOCu(atoms_obj,
                       xc = "PBE",
                       kpts = (3,3,1), 
                       encut = 400,
                       ismear = 0,
                       sigma = 0.1,
                       ediff = 1e-4,
                       prec = 'normal',
                       lcharg = False,
                       lwave = False,
                       nelmin = 4,
                       nelmdl = 6,
                       npar = 2,
                       algo = 'fast',
                       lreal = 'auto',
                       ispin = 2,
                       magmom = 0):
    """
    Takes an atoms object and sets the calculator with spin. These parameters were originally used for Jonathan Lym's Inverse Metal Oxide (ZnO/Cu(111)) runs.
    Default parameters may need to be adjusted for other systems.
    """
    formula = atoms_obj.get_chemical_formula()
    print_calc_param(calc_name = 'geometric', formula = formula, xc = xc, kpts = kpts, encut = encut, ismear = ismear, sigma = sigma, prec = prec, ediff = ediff, lcharg = lcharg,
                     lwave = lwave, nelmin = nelmin, nelmdl = nelmdl, npar = npar, algo = algo, lreal = lreal, ispin = ispin, magmom = magmom)
    calc = Vasp(xc = xc, 
            kpts = kpts, 
            encut = encut,
            ismear = ismear,
            sigma = sigma,
            ediff = ediff,
            prec = prec,
            lcharg = lcharg,
            lwave = lwave,
            nelmin = nelmin,
            nelmdl = nelmdl,
            npar = npar,
            algo = algo,
            lreal = lreal,
            ispin = ispin,
            magmom = [magmom]*len(atoms_obj))
    atoms_obj.set_calculator(calc)
    
def set_bader_calc(atoms_obj,
                   xc = "PBE",
                   kpts = (1,1,1), 
                   encut = 400,
                   ismear = 0,
                   sigma = 0.1,
                   ediff = 1e-4,
                   prec = 'normal',
                   lcharg = True,
                   lwave = False,
                   nelmin = 4,
                   nelmdl = 6,
                   npar = 2,
                   algo = 'fast',
                   lreal = 'auto',
                   ispin = 2,
                   magmom = 0,
                   laechg = True):
    formula = atoms_obj.get_chemical_formula()
    print_calc_param(calc_name = 'bader', formula = formula, xc = xc, kpts = kpts, encut = encut, ismear = ismear, sigma = sigma, prec = prec, ediff = ediff, lcharg = lcharg,
                     lwave = lwave, nelmin = nelmin, nelmdl = nelmdl, npar = npar, algo = algo, lreal = lreal, ispin = ispin, magmom = magmom)
    calc = Vasp(xc = xc, 
            kpts = kpts, 
            encut = encut,
            ismear = ismear,
            sigma = sigma,
            ediff = ediff,
            prec = prec,
            lcharg = lcharg,
            lwave = lwave,
            nelmin = nelmin,
            nelmdl = nelmdl,
            npar = npar,
            algo = algo,
            lreal = lreal,
            ispin = ispin,
            laechg = laechg,
            magmom = [magmom]*len(atoms_obj))
    atoms_obj.set_calculator(calc)
    

def set_dimer_calc_ZnOCu(atoms_obj, 
                         xc = "PBE",
                         kpts = (3,3,1),
                         encut = 400,
                         ismear = 0,
                         sigma = 0.1,
                         ediff = 1e-4,
                         lcharg = False,
                         lwave = False,
                         nelmin = 4,
                         nelmdl = 6,
                         npar = 2,
                         algo = 'fast',
                         lreal = 'auto',
                         ispin = 2,
                         magmom = 0,
                         nsw = 2000,
                         ediffg = -0.05,
                         iopt = 2,
                         ibrion = 3,
                         potim = 0,
                         ichain = 2,
                         drotmax = 6):
    """
    Takes an atoms object and sets up parameters for NEB calculations. These parameters were originally used for Jonathan Lym's Inverse Metal Oxide (ZnO/Cu(111)) runs.
    Default parameters may need to be adjusted for other systems.
    """
    formula = atoms_obj.get_chemical_formula()
    print_calc_param(calc_name = 'Dimer', formula = formula, xc = xc, kpts = kpts, encut = encut, ismear = ismear, sigma = sigma, ediff = ediff, lcharg = lcharg,
                     lwave = lwave, nelmin = nelmin, nelmdl = nelmdl, npar = npar, algo = algo, lreal = lreal, ispin = ispin, magmom = magmom)
    calc = Vasp(xc = xc, 
                kpts = kpts, 
                encut = encut,
                ismear = ismear,
                sigma = sigma,
                ediff = ediff,
                lcharg = lcharg,
                lwave = lwave,
                nelmin = nelmin,
                nelmdl = nelmdl,
                npar = npar,
                algo = algo,
                lreal = lreal,
                ispin = ispin,
                magmom = [magmom]*len(atoms_obj))

    #NEB Options
    print "-"*30
    print "NEB Calculator settings:"    
    print "nsw = %f" % nsw
    print "ediffg = %f" % ediffg
    print "iopt = %d" % iopt
    print "ibrion = %d" % ibrion
    print "potim = %d" % potim
    print "drotmax = %d" % drotmax
    print "ichain = %d" % ichain
    print "-"*30

    calc.set(nsw = nsw,
             ediffg = ediffg,
             iopt = iopt,
             ibrion = ibrion,
             potim = potim,
             drotmax = drotmax,
             ichain = ichain)
    atoms_obj.set_calculator(calc)
    return calc

def set_dimer_calc_In2O3(atoms_obj, 
                         xc = "PBE",
                         kpts = (4,3,1),
                         encut = 400,
                         ismear = 0,
                         sigma = 0.05,
                         ediff = 1e-4,
                         lcharg = False,
                         lwave = False,
                         nelmin = 4,
                         nelmdl = 6,
                         npar = 2,
                         algo = 'fast',
                         lreal = 'auto',
                         ispin = 2,
                         magmom = 0,
                         nsw = 2000,
                         ediffg = -0.05,
                         iopt = 2,
                         ibrion = 3,
                         potim = 0,
                         ichain = 2,
                         drotmax = 6):
    """
    Takes an atoms object and sets up parameters for NEB calculations. These parameters were originally used for Jonathan Lym's Inverse Metal Oxide (ZnO/Cu(111)) runs.
    Default parameters may need to be adjusted for other systems.
    """
    formula = atoms_obj.get_chemical_formula()
    print_calc_param(calc_name = 'Dimer', formula = formula, xc = xc, kpts = kpts, encut = encut, ismear = ismear, sigma = sigma, ediff = ediff, lcharg = lcharg,
                     lwave = lwave, nelmin = nelmin, nelmdl = nelmdl, npar = npar, algo = algo, lreal = lreal, ispin = ispin, magmom = magmom)
    calc = Vasp(xc = xc, 
                kpts = kpts, 
                encut = encut,
                ismear = ismear,
                sigma = sigma,
                ediff = ediff,
                lcharg = lcharg,
                lwave = lwave,
                nelmin = nelmin,
                nelmdl = nelmdl,
                npar = npar,
                algo = algo,
                lreal = lreal,
                ispin = ispin,
                magmom = [magmom]*len(atoms_obj))

    #NEB Options
    print "-"*30
    print "NEB Calculator settings:"    
    print "nsw = %f" % nsw
    print "ediffg = %f" % ediffg
    print "iopt = %d" % iopt
    print "ibrion = %d" % ibrion
    print "potim = %d" % potim
    print "drotmax = %d" % drotmax
    print "ichain = %d" % ichain
    print "-"*30

    calc.set(nsw = nsw,
             ediffg = ediffg,
             iopt = iopt,
             ibrion = ibrion,
             potim = potim,
             drotmax = drotmax,
             ichain = ichain)
    atoms_obj.set_calculator(calc)
    return calc


def set_neb_calc_ZnOCu(atoms_obj,
                       NIMAGES = 10,
                       xc = "PBE",
                       kpts = (3,3,1),
                       encut = 400,
                       ismear = 0,
                       sigma = 0.1,
                       ediff = 1e-4,
                       lcharg = False,
                       lwave = False,
                       nelmin = 4,
                       nelmdl = 6,
                       npar = 2,
                       algo = 'fast',
                       lreal = 'auto',
                       ispin = 2,
                       magmom = 0,
                       nsw = 2000,
                       ediffg = -0.10,
                       iopt = 1,
                       ibrion = 3,
                       potim = 0,
                       spring = -5,
                       lclimb = False):
    from ase.calculators.vasp_neb import Vasp as Vasp_neb
    """
    Takes an atoms object and sets up parameters for NEB calculations. These parameters were originally used for Jonathan Lym's Inverse Metal Oxide (ZnO/Cu(111)) runs.
    Default parameters may need to be adjusted for other systems.
    """
    formula = atoms_obj.get_chemical_formula()
    print_calc_param(calc_name = 'NEB', formula = formula, xc = xc, kpts = kpts, encut = encut, ismear = ismear, sigma = sigma, ediff = ediff, lcharg = lcharg,
                     lwave = lwave, nelmin = nelmin, nelmdl = nelmdl, npar = npar, algo = algo, lreal = lreal, ispin = ispin, magmom = magmom)
    calc = Vasp_neb(xc = xc, 
            kpts = kpts, 
            encut = encut,
            ismear = ismear,
            sigma = sigma,
            ediff = ediff,
            lcharg = lcharg,
            lwave = lwave,
            nelmin = nelmin,
            nelmdl = nelmdl,
            npar = npar,
            algo = algo,
            lreal = lreal,
            ispin = ispin,
            magmom = [magmom]*len(atoms_obj))


    #NEB Options
    print "-"*30
    print "Dimer Calculator settings:"    
    print "nsw = %f" % nsw
    print "ediffg = %f" % ediffg
    print "iopt = %d" % iopt
    print "ibrion = %d" % ibrion
    print "potim = %d" % potim
    print "spring = %d" % spring
    print "lclimb = %r" % lclimb
    print "-"*30

    calc.set(neb = True,
             nsw = nsw,
             ediffg = ediffg,
             iopt = iopt,
             ibrion = ibrion,
             potim = potim,
             spring = spring,
             images = NIMAGES,
             lclimb = lclimb)
    atoms_obj.set_calculator(calc)
    return calc

def set_neb_calc_In2O3(atoms_obj,
                       NIMAGES = 10,
                       xc = "PBE",
                       kpts = (4,3,1),
                       encut = 400,
                       ismear = 0,
                       sigma = 0.05,
                       ediff = 1e-4,
                       lcharg = False,
                       lwave = False,
                       nelmin = 4,
                       nelmdl = 6,
                       npar = 2,
                       algo = 'fast',
                       lreal = 'auto',
                       ispin = 2,
                       magmom = 0,
                       nsw = 2000,
                       ediffg = -0.10,
                       iopt = 1,
                       ibrion = 3,
                       potim = 0,
                       spring = -5,
                       lclimb = False):
    from ase.calculators.vasp_neb import Vasp as Vasp_neb
    """
    Takes an atoms object and sets up parameters for NEB calculations. These parameters were originally used for Jonathan Lym's Inverse Metal Oxide (ZnO/Cu(111)) runs.
    Default parameters may need to be adjusted for other systems.
    """
    formula = atoms_obj.get_chemical_formula()
    print_calc_param(calc_name = 'NEB In2O3', formula = formula, xc = xc, kpts = kpts, encut = encut, ismear = ismear, sigma = sigma, ediff = ediff, lcharg = lcharg,
                     lwave = lwave, nelmin = nelmin, nelmdl = nelmdl, npar = npar, algo = algo, lreal = lreal, ispin = ispin, magmom = magmom)
    calc = Vasp_neb(xc = xc, 
            kpts = kpts, 
            encut = encut,
            ismear = ismear,
            sigma = sigma,
            ediff = ediff,
            lcharg = lcharg,
            lwave = lwave,
            nelmin = nelmin,
            nelmdl = nelmdl,
            npar = npar,
            algo = algo,
            lreal = lreal,
            ispin = ispin,
            magmom = [magmom]*len(atoms_obj))


    #NEB Options
    print "-"*30
    print "Dimer Calculator settings:"    
    print "nsw = %f" % nsw
    print "ediffg = %f" % ediffg
    print "iopt = %d" % iopt
    print "ibrion = %d" % ibrion
    print "potim = %d" % potim
    print "spring = %d" % spring
    print "lclimb = %r" % lclimb
    print "-"*30

    calc.set(neb = True,
             nsw = nsw,
             ediffg = ediffg,
             iopt = iopt,
             ibrion = ibrion,
             potim = potim,
             spring = spring,
             images = NIMAGES,
             lclimb = lclimb)
    atoms_obj.set_calculator(calc)
    return calc


def set_calc_In2O3(atoms_obj,
                   xc = "PBE",
                   kpts = (4,3,1),
                   encut = 400,
                   ismear = 0,
                   sigma = 0.05,
                   ediff = 1e-4,
                   prec = 'normal',
                   lcharg = False,
                   lwave = False,
                   nelmin = 4,
                   nelmdl = 6,
                   npar = 2,
                   algo = 'fast',
                   lreal = 'auto',
                   ispin = 2,
                   ldautype = 2,
                   magmom = 0,
                   gamma = True):
    """
    Takes an atoms object and sets the calculator. These parameters were originally used for Jonathan Lym's In2O3 runs.
    Default parameters may need to be adjusted for other systems.
    """
    formula = atoms_obj.get_chemical_formula()
    print_calc_param(calc_name = 'relaxation for In2O3', formula = formula, xc = xc, kpts = kpts, encut = encut, ismear = ismear, sigma = sigma, prec = prec, ediff = ediff, lcharg = lcharg,
                     lwave = lwave, nelmin = nelmin, nelmdl = nelmdl, npar = npar, algo = algo, lreal = lreal, ispin = ispin, magmom = magmom, ldautype = ldautype, gamma = gamma)
    calc = Vasp(xc = xc, 
            kpts = kpts, 
            encut = encut,
            ismear = ismear,
            sigma = sigma,
            ediff = ediff,
            prec = prec,
            lcharg = lcharg,
            lwave = lwave,
            nelmin = nelmin,
            nelmdl = nelmdl,
            npar = npar,
            algo = algo,
            lreal = lreal,
            ispin = ispin,
            ldautype = ldautype,
            gamma = gamma)
    atoms_obj.set_calculator(calc)


def print_calc_param(calc_name = 'custom',
                     formula = None,
                     xc = None,
                     kpts = None, 
                     encut = None, 
                     ismear = None, 
                     sigma = None, 
                     ediff = None,
                     prec = None,
                     lcharg = None,
                     lwave = None, 
                     nelmin = None, 
                     nelmdl = None, 
                     npar = None, 
                     algo = None, 
                     lreal = None, 
                     ispin = None, 
                     magmom = None,
                     ldautype = None,
                     gamma = None):
    """
    Prints parameters used for setting up calculators.
    """
    print "-"*30
    print "Setting %s calculator to %s:" % (calc_name, formula)
    if xc is None or xc == 'PW91':
        print "Functional: PW91 [Default]"
    else:
        print "Functional: %s" % xc
        
    if kpts is None:
        print "k points: (1x1x1) [Default]"
    else:
        print "k points: (%dx%dx%d)" % kpts
    if encut is not None:
        print "Energy cutoff (eV): %f" % encut  
        
    if ismear is None or ismear == 1:
        print "Electron smearing: ismear = 1 (Method of Methfessel-Paxton order 1) [Default]"
    elif ismear == 0:
        print "Electron smearing: ismear = %d (Gaussian smearing)" % ismear
    elif ismear == -1:
        print "Electron smearing: ismear = %d (Fermi smearing)" % ismear
    elif ismear == -2:
        print "Reading from WAVECAR or INCAR file and kept fixed throughout run, ismear = %d" % ismear
    elif ismear > 0:
        print "Electron smearing: %d (Method of Methfessel-Paxton order %d)" % (ismear, ismear)
            
    if sigma is None or sigma == 0.2:
        print "Width of smearing, sigma: 0.2 [Default]"
    else:
        print "Width of smearing, Sigma: %f" % sigma

    if ediff is None or ediff == 1e-4:
        print "Convergence break condition for SC-loop, ediff: 1e-4 [Default]"
    else:
        print "Convergence break condition for SC-loop, ediff: %f" % ediff

	if prec is None or prec == 'normal':
		print "Accuracy of calculation, prec = normal [Default]" 
	else:
		print "Accuracy of calculation, prec = %s" % prec
		
    if lcharg is None or lcharg == True:
        print "LCHARG will be printed [Default]"
    else:
        print "LCHARG will NOT be printed"

    if lwave is None or lwave == True:
        print "LWAVE will be printed [Default]"
    else:
        print "LWAVE will not be printed"

    if nelmin is None or nelmin == 2:
        print "Minimum number of electronic SC steps, nelmin: 2 [Default]"
    else:
        print "Minimum number of electronic SC steps, nelmin: %d" % nelmin

    if nelmdl is None:
        print "Number of non-selfconsistent steps at the beginning, nedmdl:"
        print "Default:"
        print "-3 if ISTART=0, INIWAV = 1 and IALGO = 8"
        print "-5 if ISTART=0, INIWAV = 1 and IALGO = 48"        
    else:
        print "Number of non-selfconsistent steps at the beginning, nedmdl: %d" % nelmdl

    if npar is None:
        print "Number of bands treated in parallel, npar: Number of Cores [Default]"
    else:
        print "Number of bands treated in parallel, npar: %d" % npar

    if algo is None or algo.lower() == 'normal':    
            print "Electronic Minimization Algorithm, algo: normal (Davidson iteration scheme.) [Default]"
    elif algo.lower() == 'veryfast':
        print "Electronic Minimization Algorithm, algo: %s (RMM-DIIS iteration scheme.)" % algo
    elif algo.lower() == 'fast':
        print "Electronic Minimization Algorithm, algo: %s (Mixture of Davidson and RMM-DIIS algorithms)" % algo

    if lreal is None or lreal == False:
        print "Projection operators evaluated in reciprocal space, lreal = False [Default]"
    else:
        print "Projection operators evaluated in real space , lreal = True"
    
    if ispin is None or ispin == 1:    
        print "ispin: 1 (Spin polarization OFF) [Default]"
    else:
        print "ispin: %d (Spin polarization ON)" % ispin
        
    if magmom is None and ispin == 2:
        print "Initial magnetic momentum set to NIONS*1.0 [Default]"
    elif magmom == 0:
        print "Initial magnetic momentum set to NIONS*0.0"
    else:
        print "Initial magnetic momentum set to %r" % magmom
            
    if ldautype is None:
        print "DFT+U calculations turned OFF [Default]"
    else:
        if ldautype == 1:
            ldau_option = "Liechtenstein"
        elif ldautype == 2:
            ldau_option =  "Dudarev"
        elif ldautype == 4:
            ldau_option = "Liechtenstein (LDAU)"
        else:
            ldau_option = "Warning: Option chosen (%d) not listed" % ldautype
        print "DFT+U turned ON (%s)" % ldau_option 

    if gamma is None or gamma == False:
        print "Monkhorst-Pack scheme used to generate k points."
    else:
        print "Gamma-centered k-point sampling used to generate k points."

    print "-"*30
   
def assign_magmom(atoms_obj, ispin = None):
    #Reference: http://kitchingroup.cheme.cmu.edu/dft-book/dft.html#orgheadline8
    magmom_dict = {'H': 1.,
                   'O': 2.,
                   'Fe': 2.22,
                   'Co': 1.72,
                   'Ni': 0.61}
    magmoms = [0]*len(atoms_obj)          
    for i, atom in enumerate(atoms_obj):
        if magmom_dict.get(atom.symbol) is not None:
            magmoms[i] = magmom_dict.get(atom.symbol)
    return magmoms
