# -*- coding: utf-8 -*-
"""
Created on Thu Jun 08 07:06:52 2017

@author: Jonathan Lym
"""
from itertools import combinations, permutations
from energydiagram import ED
import numpy as np

#DFT computed energies in eV
E_slab = -394.523416
E_H2 = -6.759576
E_H2O = -14.222692
E_O2 = -9.862407


class vacancy_cluster(object):
    def __init__(self, name, energy):
        self.name = name
        self.energy = energy
        self.vacancy_sites = set(name.split('_'))
        self.n_vacancy = len(self.vacancy_sites)
                
class vacancy_clusters(object):
    def __init__(self):
        self._list = []
    
    def append(self, vacancy_cluster):
        self._list.append(vacancy_cluster)

    def extend(self, vacancy_clusters):
        self._list.extend(vacancy_clusters)

    def index(self, vacancy_sites):
        for i, vacancy_cluster in enumerate(self._list):
            if set(vacancy_sites) == set(vacancy_cluster.vacancy_sites):
                return i

    def remove(self, vacancy_sites):
        for i, vacancy_cluster in enumerate(self._list):
            if set(vacancy_sites) == set(vacancy_cluster.vacancy_sites):
                self._list.pop(i)

    def __len__(self):
        return len(self._list)
    
    def __setitem__(self, index, vacancy_cluster):
        self._list[index] = vacancy_cluster

    def __getitem__(self, index):
        return self._list[index]

    def print_vacancy_paths(self, i, ref = 'H2O', plot = False):
        #Used to determine pathway with smallest energy
        min_energy = 999.

        if ref == 'H2O':
            print 'Energies referenced to gas-phase H2O and H2.'
            E_ref = E_slab + self[i].n_vacancy * E_H2
            ref_string = '+ {}H2'.format(self[i].n_vacancy)
        elif ref == 'O2':
            print 'Energies referenced to gas-phase O2.'
            E_ref = E_slab
            ref_string = ''

        #Initializes plot
        if plot:
            diagram = ED(offset = 0.05)
            diagram.add_level(energy = 0.00, bottom_text = 'Clean {}'.format(ref_string), position = 0)
            prev_sets = []

        #Finds pathways to produce final configuration
        for permutation in permutations(self[i].vacancy_sites, self[i].n_vacancy):
            #Initial set up for the clean slab
            path = '   Clean ->'
            path_energies = 'E   0.00 '
            path_del_energies = 'ΔE       '
            
            #Energy and energy difference between stages
            E_stage = np.zeros(shape = self[i].n_vacancy+1)
            del_E = np.zeros(shape=self[i].n_vacancy)

            #Find configuration with j vacancies
            for j in range(1, self[i].n_vacancy+1):
                k = self.index(set(permutation[:j]))
                #Search for the parity if it's not found
                if k is None:
                    k = self.index(get_parity(permutation[:j]))
                    name = self[k].name
                else:
                    name = get_parity(self[k].name)

                #Record data
                path = '{} {} ->'.format(path, name)
                if ref == 'H2O':
                    E_gas = E_H2O * j + E_H2 * (self[i].n_vacancy - j)
                elif ref == 'O2':
                    E_gas = E_O2 * 0.5 * j
                E_stage[j] = self[k].energy + E_gas - E_ref
                del_E[j-1] = E_stage[j] - E_stage[j-1] 
                path_energies = '{}   {:^{}.2f}'.format(path_energies, E_stage[j], len(self[k].name))
                path_del_energies = '{}   {:^{}.2f}'.format(path_del_energies, del_E[j-1], len(self[k].name))

                if plot and (self[k].vacancy_sites not in prev_sets) and all(not np.isnan(x) for x in E_stage):
                    diagram.add_level(energy = np.round(E_stage[j], 2), position = j, bottom_text = self[k].name)
                    prev_sets.append(self[k].vacancy_sites)   

            #Print the path
            path = path[:-3] #Truncate the last arrow
            print path
            print path_energies
            print path_del_energies

            #Find the optimal path
            max_step_energy = max(del_E)
            if max_step_energy < min_energy and all([not np.isnan(E) for E in del_E]):
                min_path = path
                min_energies_string = path_energies
                min_energy = max_step_energy

        #Print most optimal path
        try:             
            print 'Min Energy Path:'
            print '#'*10
            print min_path
            print min_energies_string
            print 'Largest Energy: {:.2f} eV'.format(min_energy)
        except:
            print 'No complete paths!'
        #Plot the reaction coordinate diagram
        if plot:
            diagram.plot(unit_set = 'eV')
#            plt.savefig('{}.png'.format(self[i].name), bbox_inches = 'tight')                                     


    def print_vacancy_energy(self, site, ref = 'H2O'):
        if ref == 'H2O':
            print 'Energies referenced to gas-phase H2O and H2.'
            E_ref = E_H2O - E_H2
        elif ref == 'O2':
            print 'Energies referenced to gas-phase O2.'
            E_ref = 0.5 * E_O2
        
        print '-'*10
        print 'Vacancy formation energy of site {}'.format(site)
        print '-'*10
        print 'Vacancies Before | Vacancies After: Energy (eV)'
        for VC in self:
            if (site in VC.vacancy_sites) or (get_parity(site) in VC.vacancy_sites):
                prior_sites = VC.vacancy_sites - set([site])
                if prior_sites == VC.vacancy_sites:
                    prior_sites = VC.vacancy_sites - get_parity([site])                   
                if len(prior_sites) == 0:
                    print 'Clean: {:.2f} eV'.format(VC.energy + E_ref - E_slab)
                else:
                    i = self.index(prior_sites)
                    if i is None:
                        i = self.index(get_parity(prior_sites))
                        buf = site.replace('A', 'B')
                    else:
                        buf = site
                    prior_site_name = '_'.join(self[i].vacancy_sites)
                    print '{} | {}: {:.2f} eV'.format(prior_site_name, VC.name, VC.energy + E_ref - self[i].energy)

def get_parity(sites):
    if type(sites) is str:
        return sites.replace('A','%temp%').replace('B', 'A').replace('%temp%', 'B')
    else:                
        parity_sites = []
        for site in sites:
            if 'A' in site:
                parity_sites.append(site.replace('A', 'B'))
            elif 'B' in site:
                parity_sites.append(site.replace('B', 'A'))
        return set(parity_sites)
    
def read_vacancy_file(file_name):
    import pandas as pd

    VCs = vacancy_clusters()
    data = pd.read_excel(file_name)
    for i in range(len(data)):
        if data['Energy'][i] == 'Nan':
            energy = np.nan
        else:
            energy = data['Energy'][i]
        VCs.append(vacancy_cluster(name = str(data['Site'][i]), energy = energy))
    return VCs

def get_vacancies_from_name(string):
    return set((string.split('_')))
