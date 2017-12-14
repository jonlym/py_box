# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 11:02:54 2017

@author: Jon Lym
"""


#Files to read:
#INP.d#
#BEP
#cklink
#EAg
#EAs
#gas
#omega
#SA
#sklink
#surf
#T_flow
#thermdat
#tube
#tube_COV
#tube_mole

#OUT.d#
#Beta_out
#Ea_over_RT
#elem_comp
#EQKC_out
#gas
#general_info
#Grxn_out
#Hform_out
#Hrxn_out
#Preex_out
#rpa_vis_output
#rxn_rate_ss
#Sfrom_out
#Srxn_out
#Stoich
#surf
#tube_conv
#tube_cov_ss
#tube_gasmass_ss
#tube_gasmole_ss
#tube_mass_bal_ss
#tube_molerestart
#tube_rpa
#tube_sdot_ss

from os import chdir, system
from os.path import splitext
import subprocess
import datetime
import numpy as np
import constants as c

class chemkin_options(object):
    """
    Stores several options that will be written to tube.inp.
        reactor_type - (int or str) Reactor type:
            int    str             Description
            0      UHV/mol. beam   Ultra-high vacuum
            1      batch           Batch reactor [Default]
            2      cstr            Continuous flow stirred-tank reactor
            3      pfr             Plug flow reactor

        nruns - (int) Number of runs. Default = 1
        multiple_inputs - (boolean) If True, uses T_flow for reactor conditions. Default is True.
        lstp - (boolean) Default is True
        T - (float) Temperature of reactor in Kelvin that is used if multiple_inputs = False. Default is 298 K.
        P - (float) Pressure of reactor in atmospheres that is used if multiple_inputs = False. Default is 1 atm.
        flow_rate - (float) Volumetric flow rate in cm3/s that is used if multiple_inputs = False. Default is 1.77 cm3/s.
        SA_to_V - (float) Surface area to volume ratio in 1/cm that is used if multiple_inputs = False. Default is 6.50E2 1/cm
        t_rise - (float) Temperature rise in K used when reactor_type
        liso - (boolean) Default is True
        itpd - (int or str) 
            int    str            Description
            0      no             
            1      UHV            
            2      High Pressure  
        text - (float) Default is 923.0
        aextbyv - (float) Default is 3.571
        htc - (float) Default is 0.0
        ramp - (float) Temperature ramp in K/s.
        MARI - (string) Most abundant reaction intermediate. Symbol in thermdat
        Reactant - (string) Reactant. Symbol in thermdat
        reactor_volume - (float) Reactor volume in cm3. Default is 1 cm3
        nnodes - (int) Default is 10
        ttout - (float) Default is 0.01
        rtime - (float) Default si 100
        ntdec - (int) Default is 10
        ltra - (boolean) If set to True, transient information will be saved. If False, only steady-state information is saved. Set to False. Default is False.
        ltol - (boolean) Default is False.
        abstol - (float) Default value is 1E-10
        reltol - (float) Default value is 1E-8
        NonNeg - (boolean) If set to True, values will remain positive. Default is True.
        restart_max - (int) Maximum of times that it will be restarted. If set to <= 0, there is no limit. Default is 10
        iSolver - (boolean) If set to True, iterative solver will be used. Default is False.
        mu - (int) Upper bandwidth for Krylov solver. Default is 0
        ml - (int) Lower bandwidth for Krylov solver. Default is 0
        lcov - (boolean) Default is True.
        lStatpQ - (boolean) Default is False.
        lBEP - (boolean) Default is False.
        iScale - (int) Default is 0
        lEA - (boolean) Default is True.
        omega - (float) Default is 0.5.
        Tref_beta - (float) Default is 1
        mrpa - (int) Default is 1.
        verbose_rpa - (boolean) Default is False
        trpa - (float) Default is 900.0
        lsen - (boolean) Default is False.
        lDOE - (boolean) Default is False.
        
        catalyst_density
        site_density
    """
    def __init__(self, 
                 reactor_type = 1,
                 nruns = 1,
                 multiple_inputs = True,
                 lstp = True,
                 T = 298.,
                 P = 1.,
                 flow_rate = 1.77):
        self.nruns = nruns
        
class reaction(object):
    def __init__(self,
                 is_adsorption = False,
                 is_surface_reaction = True,
                 reaction_string = None,
                 reaction_vector = None,
                 beta = 1,
                 sticking_coeff = 0.5,
                 Ea0 = 0.,
                 Ea_over_RTs = [],
                 BEP_id = 0,
                 BEP_direction = 0,
                 thermdats = None):
        self.is_adsorption = is_adsorption
        if reaction_string is None and reaction_vector is None:
            print "Warning! No reaction string or reaction vector specified"
        else:
            if reaction_vector is not None:
                self.reaction_vector = reaction_vector
                if thermdats is not None:
                    self.reaction_string = thermdats.print_reaction(reaction_vector)
                else:
                    self.reaction_string = None
            elif reaction_string is not None:
                self.reaction_string = reaction_string
                if thermdats is not None:
                    self.reaction_vector = thermdats.reaction_to_stoich(reaction_string)
                else:
                    self.reaction_vector = None
        if is_adsorption:
            self.sticking_coeff = sticking_coeff
        self.is_surface_reaction = is_surface_reaction
        self.beta = beta
        self.sticking_coeff = sticking_coeff
        self.Ea0 = Ea0
        self.Ea_over_RTs = Ea_over_RTs
        self.BEP_id = BEP_id
        self.BEP_direction = BEP_direction                
        
class reactions(object):
    
    def __init__(self, reactions = []):
        self._reactions = []
        for reaction in reactions:
            self.append(reaction)
            
    def append(self, reaction):
        self._reactions.append(reaction)

    def extend(self, reactions):
        self._reactions.extend(reactions)

    def remove(self, index):
        self._reactions.pop(index)

    def __len__(self):
        return len(self._reactions)
    
    def __setitem__(self, index, reaction):
        self._reactions[index] = reaction

    def __getitem__(self, index):
        return self._reactions[index]
        
    def __str__(self):
        rxn = ''
        for reaction in self:
            if reaction.reaction_string is not None:
                rxn += "%r\n" % reaction.reaction_string
        return rxn

class BEP(object):
    def __init__(self,
                 BEP_type = 0.,
                 m = 0.,
                 b = 0.,
                 BEP_direction = -1.,
                 mean_uncertainty = 0.,
                 std_uncertainty = 0.,
                 description = ''):
        self.BEP_type = BEP_type
        self.m = m
        self.b = b
        self.BEP_direction = BEP_direction
        self.mean_uncertainty = mean_uncertainty
        self.std_uncertainty = std_uncertainty
        self.description = description


class BEPs(object):
    def __init__(self, BEPs = []):
        self._BEPs = []
        for BEP in BEPs:
            self.append(BEP)
            
    def append(self, BEP):
        self._BEPs.append(BEP)

    def extend(self, BEPs):
        self._BEPs.extend(BEPs)

    def remove(self, index):
        self._BEPs.pop(index)

    def __len__(self):
        return len(self._BEPs)
    
    def __setitem__(self, index, BEP):
        self._BEPs[index] = BEP

    def __getitem__(self, index):
        return self._BEPs[index]
        
    def print_BEPs(self):
        """
        Prints a summary of the thermdat list.
        """
        print "Index\tSlope\tIntercept (kcal/mol)\tDescription"
        for i, BEP in enumerate(self):
            print "%d\t%f\t%f\t%s" % (i, BEP.m, BEP.b, BEP.description)



def read_reactions(in_path):
    reactions_list = reactions()
    with open(in_path, 'r') as reaction_file:
        for line in reaction_file:
            if line[0] != '!':
                #Split data and remove unnecessary characters
                data = line.split(',')
                reaction_string = data[0]
                if data[1] == '' or data[1] == '\n':
                    is_adsorption = False
                else:
                    if 'true' in data[1].lower():
                        is_adsorption = True
                    elif 'false' in data[1].lower():
                        is_adsorption = False

                if data[2] == '' or data[2] == '\n':
                    is_surface_reaction = True
                else:
                    if 'true' in data[2].lower():
                        is_surface_reaction = True
                    elif 'false' in data[2].lower():
                        is_surface_reaction = False

                if data[3] == '' or data[3] == '\n':
                    beta = 1
                else:
                    beta = float(data[3])

                if data[4] == '' or data[4] == '\n':
                    sticking_coeff = 0.5
                else:
                    sticking_coeff = float(data[4])

                if data[5] == '' or data[5] == '\n':
                    BEP_id = 0
                else:                
                    BEP_id = int(float(data[5]))
                    
                if data[6] == '' or data[6] == '\n':
                    BEP_direction = 0
                else:
                    BEP_direction = int(float(data[6]))

                if data[7] == '' or data[7] == '\n':
                    Ea0 = 0
                else:                
                    Ea0 = float(data[7])

                Ea_over_RTs = [float(i) for i in data[8:] if i != '' and i != '\n']
                reactions_list.append(reaction(is_adsorption = is_adsorption,
                                               is_surface_reaction = is_surface_reaction,
                                               reaction_string = reaction_string,
                                               beta = beta,
                                               sticking_coeff = sticking_coeff,
                                               Ea0 = Ea0,
                                               Ea_over_RTs = Ea_over_RTs,
                                               BEP_id = BEP_id,
                                               BEP_direction = BEP_direction))
    return reactions_list

def read_BEPs(in_path):
    BEPs_list = BEPs()
    with open(in_path, 'r') as BEP_file:
        for line in BEP_file:
            if line[0] != '!':
                #Split data and remove unnecessary characters
                data = line.split(',')
                BEP_type = int(float(data[0]))

                if data[1] == '' or data[1] == '\n':
                    m = 0.
                else:
                    m = float(data[1])

                if data[2] == '' or data[2] == '\n':
                    b = 0.
                else:
                    b = float(data[2])

                if data[3] == '' or data[3] == '\n':
                    BEP_direction = 0
                else:
                    BEP_direction = int(float(data[3]))

                if data[4] == '' or data[4] == '\n':
                    mean_uncertainty = 0.
                else:                
                    mean_uncertainty = float(data[4])
                    
                if data[5] == '' or data[5] == '\n':
                    std_uncertainty = 0.
                else:
                    std_uncertainty = float(data[5])

                if data[6] == '' or data[6] == '\n':
                    description = ''
                else:
                    description = data[6].replace('\n', '')

                BEPs_list.append(BEP(BEP_type = BEP_type, 
                                    m = m, 
                                    b = b, 
                                    BEP_direction = BEP_direction, 
                                    mean_uncertainty = mean_uncertainty,
                                    std_uncertainty = std_uncertainty,
                                    description = description))
    return BEPs_list


def extract_chemkin(path):
    """
    Main script that will read files in INP.d and OUT.d and export the values to a spreadsheet.
    """
    chdir(path)

def run_chemkin(log_file = 'output.log', clean_output = True):
    script_path = './chemkin.sh'
    p = subprocess.Popen([script_path], stdin = subprocess.PIPE)
    p.communicate(input='Y\nY\nN\n')
    clean_terminal_output(in_path = log_file)
    
def clean_terminal_output(in_path, out_path = None, print_line_count = True):
    """
    If the terminal output was saved (e.g. by using the nohup command) then this script will summarize lines unnecessarily
    increase the file size.
    """

    i = 0
    del_line_i = False
    del_line_j = False
    #Terms in this list will be summarized
    clean_terms = ['DASPK--  AT CURRENT T (=R1)', 
                   'DASPK--  TAKEN ON THIS CALL BEFORE REACHING TOUT',
                   'In above message,  R1 =']

    if out_path == None:
        in_file, in_ext = splitext(in_path)
        out_path = '%s%s%s' % (in_file, '_summary', in_ext)
        
    print "Reading from %s and writing to %s" % (in_path, out_path)
    input_file = open(in_path, 'r')
    output_file = open(out_path, 'w')
    for line in input_file:
        #Determine if the line should be deleted
        for clean_term in clean_terms:
            if clean_term in line:
                del_line_i = True
                i += 1
                break
        #If the current line does not need to be deleted
        if not del_line_i:
            #But the previous line does
            if del_line_j:
                #Write a summary of how many lines were deleted
                if print_line_count:
                    output_file.write('Truncated %d lines\n' % i)
            #Then write the current line
            output_file.write(line)
            i = 0
        #Reset the booleans
        del_line_j = del_line_i
        del_line_i = False            
    input_file.close()
    output_file.close()
    print "Completed chemkin.clean_terminal_output"

def write_gas(thermdats, reactions, out_path = 'gas.inp', metal_site = 'CU'):
    with open(out_path, 'w') as gas:
        element_header = """!File generated on %s using chemkin.write_gas() by Jonathan Lym
ELEMENTS
O
H
C
N
%s
END

SPECIES
""" % (datetime.datetime.now().strftime("%m/%d/%Y %H:%m"), metal_site)
        gas.write(element_header)
        for thermdat in thermdats:
            if thermdat.is_gas:
                gas.write('%s\n' % thermdat.symbol)
        reaction_header = """END
!THERMO
! Insert GRI-Mech thermodynamics here or use in default file
!END
REACTIONS
"""
        gas.write(reaction_header)
        for reaction in reactions:
            if not reaction.is_surface_reaction:
                reaction_string = reaction.reaction_string.replace('=', '<=>').replace(' ', '')
                gas.write('%s%s%.3E    %.3f    %.2f\n' % (reaction_string, ' '*(41 - len(reaction_string)), c.kb('J/K')/c.h('J s'), reaction.beta, c.convert_unit(from_ = 'kcal', to = 'cal')*reaction.Ea0))
        gas.write('END')
    
def write_surf(thermdats, reactions, out_path = 'surf.inp', site_density = 2.98897E-9, metal_density = 21.4, metal_site = 'CU'):
    with open(out_path, 'w') as surf:
        species_header = """!File generated on %s using chemkin.write_surf() by Jonathan Lym
!Surface species
!==============================================================================
!Format is Label(S)/number of sites occupied/

SITE/SURFACE/      SDEN/%.5E/

""" % (datetime.datetime.now().strftime("%m/%d/%Y %H:%m"), site_density)
        surf.write(species_header)
        for thermdat in thermdats:
            if not thermdat.is_gas:
                #Ignore the metal species until later
                if np.sum(thermdat.CHON) != 0:
                    surf.write('  %s/%d/\n' % (thermdat.symbol, thermdat.site_type))
                else:
                    if '(S)' in thermdat.symbol:
                        empty_site = thermdat.symbol
                    elif '(B)' in thermdat.symbol:
                        bulk_site = thermdat.symbol
                    else:
                        print "Warning: %s does not contain the number of C, H, O, or N" % thermdat.symbol
        surf.write('  %s\n\nBULK %s/%.1f/\nEND\n' % (empty_site, bulk_site, metal_density))
        
        reaction_header = """!Reactions
!==============================================================================

REACTIONS MWOFF KCAL/MOLE

!***************  Notes about the format of this surf.inp file. ***************
!Entries have the format <Balanced Eq.> <Pre-exp./Stick. Coeff.> <Beta> <EA>.
!Modified Arrhenius kinetics are assumed: k=A*((T/To)**beta)*exp(-EA/RT).
!In this surf.inp, To=1 K and A=(kB/h)*(sden)**(1-m) where m is the number of
!adsorbed reactants in the balanced equation.
"""
        surf.write(reaction_header)                            
        for reaction in reactions:
            if reaction.is_surface_reaction:
                if reaction.reaction_string is None:
                    reaction.reaction_string = thermdats.print_reaction(reaction.reaction_vector)
                if reaction.is_adsorption:
                    A = reaction.sticking_coeff
                    stick = '\n  STICK'
                else:
                    if reaction.reaction_vector is None:
                        reaction.reaction_vector = thermdats.reaction_to_stoich(reaction.reaction_string)
                    num_adsorbed_species = 0
                    for i, coefficient in enumerate(reaction.reaction_vector):
                        #If this is a surface species (includes empty sites)
                        if coefficient < 0 and not thermdats[i].is_gas and '(B)' not in thermdats[i].symbol:
                            num_adsorbed_species -= coefficient
                    A = c.kb('J/K')/(c.h('J s') * site_density ** (num_adsorbed_species-1))
                    stick = ''
                surf.write('  %s%s%.2E %.2f %.2f%s\n' % (reaction.reaction_string, ' '*(59 - len(reaction.reaction_string)), A, reaction.beta, reaction.Ea0, stick))
        surf.write('END')
        
def write_BEP(reactions, BEPs, out_path = 'BEP.inp'):
    BEP_header = """!File generated on %s using chemkin.write_BEP() by Jonathan Lym
!BEP input.  In the lines below specify the number and specifics of linear correlations used.
!------------------------------------------------
   %d          !number of correlations
!------------------------------------------------
! BEP #, BEP Type, m, b, decomposition/synthesis, mean/std dev uncertainty
!   BEP Type: -1/0/1 is TSS IS/BEP/TSS FS
!   m/b defined by
!     -1: (ETS)=m(EIS)+b 
!      0: (Ea)=m(deltaHrxn)+b
!      1: (ETS)=m(EFS)+b 
!     All intercepts are in kcal/mol and EIS, EFS, ETS are the heats of
!       formation of the species (/NOT/ the binding energies)
!   decomposition/synthesis reference direction denoted by -1/1
!   If uncertainty values not needed, can use zeros as dummy entries
!------------------------------------------------
""" % (datetime.datetime.now().strftime("%m/%d/%Y %H:%m"), len(BEPs))
    with open(out_path, 'w') as BEP_file:
        BEP_file.write(BEP_header)
        for i, BEP in enumerate(BEPs):
            if BEP.BEP_direction < 0:
                space = ''
            else:
                space = ' '
            BEP_file.write('  %d  %d  %.4f  %.2f %s%d  %.1f  %.1f    !%s\n' % (i+1, BEP.BEP_type, BEP.m, BEP.b, space, BEP.BEP_direction, BEP.mean_uncertainty, BEP.std_uncertainty, BEP.description))
        
        reaction_count = 0
        for reaction in reactions:
            if reaction.BEP_id != 0:
                reaction_count += 1
        reaction_header = """!------------------------------------------------
!Specify BEP number and reaction direction for all reactions below.
!The BEP number is given above (enter 0 if BEP should not be used for
!this reaction). The reaction direction is one of -1/1 for
!decomposition/synthesis. If BEPs are not used, then the reaction
!direction does not need to be specified (use 0 as a dummy value).
!Only those reactions using BEPs need to be specified. This is done
!via the reaction string in the third column. This reaction string must
!/exactly/ match the reaction string in surf.out (aside from the line break
!and writing repeated species with a summed stoichiometric coefficient);
!any deviation is a fatal error which will terminate the program. The first
!line following the comments should be the total number of reactions using
!BEP correlations.
!------------------------------------------------
%d   Number of non-zero values
""" % reaction_count
        BEP_file.write(reaction_header)
        for reaction in reactions:
            reaction_string = reaction.reaction_string.replace('=', '<=>').replace(' ', '')
            if reaction.BEP_id != 0:
                if reaction.BEP_direction < 0:
                    space = ''
                else:
                    space = ' '
                BEP_file.write("%d %s%d     %s\n" % (reaction.BEP_id, space, reaction.BEP_direction, reaction_string))
        BEP_file.write('EOF')