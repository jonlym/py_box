import py_box.cluster_expansion as CE
from random import shuffle
from copy import copy

def get_bins(data, k = None, shuffle_list = True):
    if k is None:
        k = len(data)

    list_copy = copy(data)
    if shuffle_list:
        shuffle(list_copy)
    n_element = len(data)/k

    if type(data) is CE.Configurations:
        bins = CE.Configurations()
    elif type(data) is CE.Clusters:
        bins = CE.Clusters()
    elif type(data) is list:
        bins = []
    else:
        print 'Bin data is {}'.format(type(data))
    bins.extend([list_copy[i:i + n_element] for i in xrange(0, len(data), n_element)])

    #If more bins than requested, distribute elements in the last bin
    while len(bins) > k:
        for j, element in enumerate(bins[-1]):
            bins[j%k].append(element)
        bins.pop(-1)
    return bins

def separate_bins(i, data):
    # test_set = copy(data[i])
    #
    if type(data) is CE.Configurations:
        training_set = CE.Configurations()
        test_set = CE.Configurations()
    elif type(data) is CE.Clusters:
        training_set = CE.Clusters()
        test_set = CE.Clusters()
    elif type(data) is list:
        training_set = []
        test_set = []
    test_set.extend(copy(data[i]))
    for j in xrange(len(data)):
        if j != i:
            training_set.extend(copy(data[j]))
    return (training_set, test_set)
