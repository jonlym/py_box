from py_box3 import parse_formula
import numpy as np

def get_query_criteria_from_elements(elements):
    """
    Takes an elements list and converts it to the pymatgen criteria

    Arguments
        elements - list
            List of elements
            e.g. ['Al', 'O']
    Returns
        str
            pymatgen criteria
            e.g. 'Al-O'
    """
    return '-'.join(elements)

def get_unit_energy(query_result):
    """
    Gets the energy per 'pretty formula' unit
    
    Parameters
        query_result - dict
            Dictionary from pymatgen query that contains the following fields:
                - pretty_formula
                - unit_cell_formula
                - energy
    Returns
        float
            Unit energy per 'pretty formula' unit
    """
    pretty_formula = parse_formula(query_result['pretty_formula'])
    #Calculate number of units in cell
    n_units  = np.zeros(len(query_result['unit_cell_formula']))
    for i, (element, unit_val) in enumerate(query_result['unit_cell_formula'].items()):
        n_units[i] = unit_val/pretty_formula[element]

    #Warn if the number of units is not equal for all elements
    if not np.all([n_units == n_units[0]]):
        raise ValueError('Number of units not consistent for all elements.')

    unit_oxide_energy = query_result['energy'] / n_units[0]
    return unit_oxide_energy
