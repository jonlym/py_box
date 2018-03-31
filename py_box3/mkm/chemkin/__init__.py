# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 14:57:39 2016

@author: Jonathan Lym
"""

from py_box3.constants import T0, convert_unit
from ase.io import read
from py_box3.thermo.thermdat import Thermdat
from py_box3.thermo.thermdats import Thermdats
import numpy as np

class Chemkin(object):
    def __init__(self, 
                species = None,
                sites = None,
                reactions = None,
                BEPs = None,
                LSRs = None, 
                DOEs = None,
                GAs = None,
                SAs = None,
                StatpQ = None,
                reactor_type = 1,
                n_runs = 1,
                multi_input = True,
                standard_T_and_P = True,
                Ts = [],
                Ps = [],
                Qs = [],
                SA_Vs = [],
                T_rise = 0.,
                isothermal = True,
                linear_T_ramp = False,
                external_T = 923.,
                heat_transfer_area_to_volume = 3.571,
                heat_transfer_coefficient = 0.,
                TPD_ramp = 0.,
                MARI = '',
                reactant = '',
                volume = 100.,
                nnodes = 1,
                ttout = 1.e-2,
                rtime = 1.e-4,
                ntdec = 10.,
                save_transient = False,
                set_equation_tolerance = True,
                absolute_tolerance = 1.e-10,
                relative_tolerance = 1.e-8,
                non_negative_composition = True,
                restart_max = 0,
                use_iterative_solver = False,
                upper_bandwidth = 0,
                lower_bandwidth = 0,
                use_coverage_effects = False,
                use_binding_energy_corrections = False,
                use_BEPs = False,
                use_LSRs = False,
                use_different_activation_energy = False,
                use_omega = False,
                omega = 0.5,
                T_ref = 1.,
                reaction_path_analysis_mode = 1,
                verbose_reaction_path_analysis = False,
                reaction_path_analysis_T = 900.,
                sensitivity_analysis = False,
                design_of_experiments = False):
        #Objects
        self.species = species
        self.sites = sites
        self.reactions = reactions
        self.BEPs = BEPs
        self.LSRs = LSRs
        self.DOEs = DOEs
        self.GAs = GAs
        self.SAs = SAs
        self.StatpQ = StatpQ

        #Reactor parameters
        self.reactor_type = reactor_type
        self.n_runs = n_runs
        self.multi_input = multi_input
        self.standard_T_and_P = standard_T_and_P
        self.Ts = Ts
        self.Ps = Ps
        self.Qs = Qs
        self.SA_Vs = SA_Vs
        self.T_rise = T_rise
        self.external_T = external_T
        self.heat_transfer_area_to_volume = heat_transfer_area_to_volume
        self.heat_transfer_coefficient = heat_transfer_coefficient
        self.TPD_ramp = TPD_ramp
        self.MARI = MARI
        self.reactant = reactant
        self.volume = volume

        #Reactor Options
        self.isothermal = isothermal
        self.linear_T_ramp = linear_T_ramp

        #Solver options
        self.nnodes = nnodes
        self.ttout = ttout
        self.rtime = rtime
        self.ntdec = ntdec
        self.save_transient = save_transient
        self.set_equation_tolerance = set_equation_tolerance
        self.absolute_tolerance = absolute_tolerance
        self.relative_tolerance = relative_tolerance
        self.non_negative_composition = non_negative_composition
        self.restart_max = restart_max
        self.use_iterative_solver = use_iterative_solver
        self.upper_bandwidth = upper_bandwidth
        self.lower_bandwidth = lower_bandwidth

        #Reaction options
        self.use_coverage_effects = use_coverage_effects
        self.use_binding_energy_corrections = use_binding_energy_corrections
        self.use_BEPs = use_BEPs
        self.use_LSRs = use_LSRs
        self.use_different_activation_energy = use_different_activation_energy
        self.use_omega = use_omega
        self.omega = omega
        self.T_ref = T_ref

        #Output options
        self.reaction_path_analysis_mode = reaction_path_analysis_mode
        self.verbose_reaction_path_analysis = verbose_reaction_path_analysis
        self.reaction_path_analysis_T = reaction_path_analysis_T
        self.sensitivity_analysis = sensitivity_analysis
        self.design_of_experiments = design_of_experiments

    @classmethod
    def from_INP(self, path = '.'):

        sites = Sites.from_surf_inp(path = os.path.join(path, 'surf.inp'))
        species = Species.from_thermdat(path = os.path.join(path, 'thermdat'))
        species.get_sites(path = os.path.join(path, 'surf.inp'))

        gas_reactions = Reactions.from_gas_inp(path = os.path.join(path, 'gas.inp'))
        surf_reactions = Reactions.from_surf_inp(path = os.path.join(path, 'surf.inp'))
        reactions = copy(gas_reactions).extend(copy(surf_reactions))

        input_dict = self.read_tube_inp(path = os.path.join(path, 'tube.inp'), return_dict = True)

        #Optional Objects
        if tube_dict['use_BEPs']:
            input_dict['BEPs'] = BEPs.from_BEP_inp(path = os.path.join(path, 'BEP.inp'))

        if tube_dict['use_LSRs']:
            input_dict['LSRs'] = LSRs.from_Scale_inp(path = os.path.join(path, 'Scale.inp'))

        if tube_dict['design_of_experiments']:
            input_dict['DOEs'] = DOEs.from_DOE_inp(path = os.path.join(path, 'DOE.inp'))

        if tube_dict['use_GAs']:
            input_dict['GAs'] = GAs.from_GA_inp(path = os.path.join(path, 'GA.inp'))

        if tube_dict['sensitivity_analysis']:
            input_dict['SAs'] = SAs.from_SA_inp(path = os.path.join(path, 'SA.inp'))

        if tube_dict['use_binding_energy_corrections']:
            input_dict['StatpQ'] = StatpQ.from_StatpQ_inp(path = os.path.join(path, 'StatpQ.inp'))

        if tube_dict['multi_input']:
            (Ts, Ps, Qs, SA_Vs) = self.read_T_flow_inp(path = os.path.join(path, 'T_flow.inp'))
            if tube_dict['use_different_activation_energy']:
                reactions.read_EAs_inp(path = os.path.join(path, 'EAs.inp'))
                reactions.read_EAg_inp(path = os.path.join(path, 'EAg.inp'))

        return cls(species = species, sites = sites, reactions = reactions, **input_dict)

    def read_tube_inp(self, path = 'tube.inp', return_dict = True):
        tube_dict = dict()
        with open(path, 'r') as f_ptr:
            i = 0
            for line in f_ptr:
                #Skip lines
                if '!' == line[0] or 'EOF' in line:
                    continue

                data = [x for x in line.replace('\n', '').split(' ') if x != '']
                if i == 0:
                    tube_dict['reactor_type'] = int(data[0])
                    tube_dict['n_runs'] = int(data[1])
                    tube_dict['multi_input'] = char_to_boolean(data[2])
                elif i == 1:
                    tube_dict['standard_T_and_P'] = char_to_boolean(data[0])
                    tube_dict['Ts'] = [float(data[1])]
                    tube_dict['Ps'] = [float(data[2])]
                    tube_dict['Qs'] = [float(data[3])]
                    tube_dict['SA_Vs'] = [float(data[4])]
                    tube_dict['T_rise'] = float(data[5])
                elif i == 2:
                    tube_dict['isothermal'] = char_to_boolean(data[0])
                    tube_dict['linear_T_ramp'] = int(data[1])
                elif i == 3:
                    tube_dict['external_T'] = float(data[0])
                    tube_dict['heat_transfer_area_to_volume'] = float(data[1])
                    tube_dict['heat_transfer_coefficient'] = float(data[2])
                    tube_dict['TPD_ramp'] = float(data[3])
                elif i == 4:
                    tube_dict['MARI'] = data[0]
                    tube_dict['reactant'] = data[1]
                elif i == 5:
                    tube_dict['volume'] = float(data[0])
                    tube_dict['nnodes'] = int(data[1])
                    tube_dict['ttout'] = float(data[2])
                    tube_dict['rtime'] = float(data[3])
                    tube_dict['ntdec'] = int(data[4])
                    tube_dict['save_transient'] = char_to_boolean(data[5])
                elif i == 6:
                    tube_dict['set_equation_tolerance'] = char_to_boolean(data[0])
                    tube_dict['absolute_tolerance'] = float(data[1])
                    tube_dict['relative_tolerance'] = float(data[2])
                    tube_dict['non_negative_composition'] = char_to_boolean(data[3])
                    tube_dict['restart_max'] = int(data[4])
                elif i == 7:
                    if data[0] == '0':
                        tube_dict['use_iterative_solver'] = False
                    elif data[0] == '1':
                        tube_dict['use_iterative_solver'] = True
                    else:
                        raise Exception('Invalid value for iSolver, {}'.format(data[0]))
                    tube_dict['upper_bandwidth'] = int(data[1])
                    tube_dict['lower_bandwidth'] = int(data[2])
                elif i == 8:
                    tube_dict['use_coverage_effects'] = char_to_boolean(data[0])
                    tube_dict['use_binding_energy_corrections'] = char_to_boolean(data[1])
                    tube_dict['use_BEPs'] = char_to_boolean(data[2])
                    if data[3] == '0':
                        tube_dict['use_LSRs'] = False
                    elif data[3] == '3':
                        tube_dict['use_LSRs'] = True
                    else:
                        raise Exception('Invalid value for iScale, {}'.format(data[3]))
                    tube_dict['use_different_activation_energy'] = char_to_boolean(data[4])
                    tube_dict['use_omega'] = char_to_boolean(data[5])
                    tube_dict['omega'] = float(data[6])
                    tube_dict['T_ref'] = float(data[7])
                elif i == 9:
                    tube_dict['reaction_path_analysis_mode'] = int(data[0])
                    tube_dict['verbose_reaction_path_analysis'] = char_to_boolean(data[1])
                    tube_dict['reaction_path_analysis_T'] = float(data[2])
                    tube_dict['sensitivity_analysis'] = char_to_boolean(data[3])
                    tube_dict['design_of_experiments'] = char_to_boolean(data[4])
                i += 1
        return tube_dict

    def write_tube_inp(self, path = 'tube.inp'):
        lines = []
        lines.append('!irxtr (0=UHV/mol. beam, 1=batch, 2=cstr, 3=pfr)    nruns  MultiInput')
        #lines.append('{}{}{}{}{}'.format(self.reactor_type))
        lines.append('!lstp  t[K]   p[atm]  velo[cm3/s]  abyv[cm-1]  trise[K]')
        lines.append('!liso(yes=T,no=F) itpd (0=no, 1=UHV, 2=High Pressure) (itpd overrides liso)')
        lines.append('!text   aextbyv htc  ramp [K/s]')
        lines.append('!MARI               Reactant')
        lines.append('!rlen[cm3]  nnodes ttout [s] rtime [s]  ntdec  ltra (F=only SS saved, T=transient saved)')
        lines.append('!ltol  abstol reltol  NonNeg(F/T: constraints off/on) restart_max (<=0 means no limit)')
        lines.append('!iSolver (0/1: iterative solver off/on)  mu  ml (upper/lower bandwidths for Krylov solver)')
        lines.append('!lcov lStatpQ lBEP iScale lEA lomega omega Tref_beta (0: Tref=300K; 1: Tref=1K)')
        lines.append('!mrpa verbose_rpa trpa    lsen   lDOE')
        lines.append('EOF')
        with open(path, 'w') as f_ptr:
            f_ptr.write(lines[0])




def char_to_boolean(character):
    if character.lower() == 't':
        return True
    elif character.lower() == 'f':
        return False
    else:
        raise Exception('Invalid character, {}'.format(character))

def boolean_to_char(boolean):
    if boolean:
        return 'T'
    else:
        return 'F'

def read_freq(freq_file_path, freq_cut_off = 0, verbose = True, warn = True):
    """
    Reads the csv file that has the frequencies of the species.
    The CSV file must be of the format:
    Symbol, is_gas, #C, #H, #O, #N, H0, potentialenergy, geometry, symmetry, frequencies
    """
    if verbose:
        print(("Reading from file: %s" % freq_file_path))
    thermdats = Thermdats()
    freq_file = open(freq_file_path, 'r')
    with open(freq_file_path, 'r') as freq_file:
        for line in freq_file:
            if line[0] is not '!':
                #Split data and remove unnecessary characters
                data = line.split(',')
                #Sorting data into respective sets
                #Symbol
                symbol = data[0]
                #is_gas            
                is_gas = int(data[1])            
                #CHON            
                CHON = [int(float(i)) for i in data[2:6]]
                if data[6] == '':
                    data[6] = '0'
                H0 = float(data[6])
                #potentialenergy            
                if data[7] == '':
                    data[7] = '0'
                potentialenergy = float(data[7])            
                #geometry
                if data[8] == '':
                    geometry = None
                else:
                    geometry = data[8]
                #symmetrynumber
                if data[9] == '':
                    symmetrynumber = None
                else:
                    symmetrynumber = int(data[9])
                #spin
                if data[10] == '':
                    spin = None
                else:
                    spin = float(data[10])
                #atoms
                if data[11] == '':
                    atoms = None
                else:
                    atoms = read(data[11])
                #vib_freq
                vib_freq = np.array([float(i) for i in data[12:] if i != '' and i != '\n'])
                if verbose:
                    print(("Ignoring frequencies below %f cm^-1" % freq_cut_off))
                vib_freq = vib_freq[vib_freq >= freq_cut_off]
                if verbose:
                    print(("Importing %s" % symbol))
                thermdats.append(Thermdat(symbol = symbol,
                                              is_gas = is_gas,
                                              CHON = CHON,
                                              H0 = H0,
                                              potentialenergy = potentialenergy,
                                              geometry = geometry,
                                              symmetrynumber = symmetrynumber,
                                              spin = spin,
                                              atoms = atoms,                                            
                                              vib_freq = vib_freq,
                                              verbose = verbose,
                                              warn = warn))
    return thermdats                               

def read_ref(path, verbose = True, warn = True):
    #Read the reference files
    species_ref = read_freq(path, verbose = verbose, warn = warn)
    
    #Prepare the fundamental species
    e_list = ['C', 'H', 'O', 'N']
    rm_list = []
    fund_species_CHON = np.array([[1, 0, 0, 0],
                                 [0, 2, 0, 0],
                                 [0, 0, 2, 0],
                                 [0, 0, 0, 2]])
    
    n_s = len(species_ref)
    H0_dft = np.array([0.]*n_s)
    H0_exp = np.array([0.]*n_s)
    species_CHON = np.array([[0]*4]*n_s)
    e_sum_list = np.array([0.]*4)
    
    for i, species in enumerate(species_ref):
        #Calculate properties at T0
        H0_dft[i] = species.IdealGasThermo.get_enthalpy(T0('K'), verbose = False)
        H0_exp[i] = species.H0*convert_unit(from_ ='kJ/mol', to = 'eV/molecule') #species.H0 from literature

        #Sets up the CHON matrix
        species_CHON[i, :] = species.CHON
        e_sum_list += species.CHON #Used to determine if elements are not needed. e.g. O, N not needed for hydrocarbons    

    #Cleans up the CHON matrix if necessary    
    for i, (element, e_sum) in reversed(list(enumerate(zip(e_list, e_sum_list)))):
        if e_sum == 0:
            if verbose:
                print(('Element %s not needed' % element))
            species_CHON = np.delete(species_CHON, i, 1)
            fund_species_CHON = np.delete(fund_species_CHON, i, 0)
            fund_species_CHON = np.delete(fund_species_CHON, i, 1)
            rm_list.append(i)
    n_e = len(fund_species_CHON)
    if n_e != n_s and verbose:
        print(("Warning: Mismatch between number of elements (%d) and number of reference species." % (n_e, n_s)))
                    
    #Finds the fundamental energy
    H0_fund = np.array([0.]*n_e) #Formation enthalpy for fundamental species. Usually 0
    ref_species_from_fund = np.dot(species_CHON, np.linalg.inv(fund_species_CHON)) #Matrix to go from fundamental species --> reference species
    H0_rxn = H0_exp - np.dot(ref_species_from_fund, H0_fund) #Since H0_fund = 0, this is H0_exp
    energies_fund = np.linalg.solve(ref_species_from_fund, (H0_dft-H0_rxn)) #Energy of fundamental species. Used as correction factor
    return (energies_fund, fund_species_CHON, rm_list)