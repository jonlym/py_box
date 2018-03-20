import numpy as np

index_dict = {'O1A': 0,
              'O2A': 1,
              'O3A': 2,
              'O4A': 3,
              'O5A': 4,
              'O6A': 5,
              'O1B': 6,
              'O2B': 7,
              'O3B': 8,
              'O4B': 9,
              'O5B': 10,
              'O6B': 11}

pos_dict = {0: 75,
            1: 76,
            2: 46,
            3: 73,
            4: 43,
            5: 68,
            6: 74,
            7: 77,
            8: 47,
            9: 72,
            10: 42,
            11: 69}

def get_interactions_from_sites(sites, delimiter = '_'):
    if type(sites) is str:
        sites = sites.split(delimiter)
    buf1 = index_dict[sites[0]]
    buf2 = index_dict[get_parity(sites[0])]
    for site in sites[1:]:
        buf1 = '{}, {}'.format(buf1, index_dict[site])
        buf2 = '{}, {}'.format(buf2, index_dict[get_parity(site)])
    # #Take out ending ', '
    # buf1 = buf1[:-2]
    # buf2 = buf2[:-2]

    return '{} | {}'.format(buf1, buf2)

def get_sigma_from_sites(sites, delimiter = '_'):
    sigma = np.ones(shape = (len(index_dict)))

    if sites.lower() == 'clean':
        return sigma

    if type(sites) is str:
        sites = sites.split(delimiter)
    for site in sites:
        sigma[index_dict[site]] = -1
    return sigma

def get_vacancy_indices_from_sigma(sigma):
  return [pos_dict[i] for i, x in enumerate(sigma) if x == -1]


def get_parity(site):
    #If site is a name
    if type(site) is str:
        return site.replace('A', '%temp%').replace('B', 'A').replace('%temp%', 'B')
    #If site is an index
    elif type(site) is int:
        if site > 5:
            return site-6
        else:
            return site+6

