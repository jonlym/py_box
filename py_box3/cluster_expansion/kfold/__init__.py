import py_box3.cluster_expansion as CE
from random import shuffle
from copy import copy

def get_bins(data, k = None, shuffle_data = True):
    data_copy = copy(data)

    if k is None:
        k = len(data)
    elif shuffle_data:
        shuffle(data_copy)

    n_element = len(data)/k

    bins = _return_empty_structure(data)
    bins.extend([data_copy[i:i + n_element] for i in range(0, len(data), n_element)])

    #If more bins than requested, distribute elements in the last bin
    while len(bins) > k:
        for j, element in enumerate(bins[-1]):
            bins[j%k].append(element)
        bins.remove(index=-1)
    return bins

def separate_bins(i, data):
    training_set = _return_empty_structure(data)
    test_set = _return_empty_structure(data)

    # if all([len(x) == 1 for x in data]):
    #     #If bins only contain one element, use leave-out-one method
    #     test_set.extend(copy(data))
    # else:
    #     #Otherwise, do k fold cross validation
    test_set.extend(copy(data[i]))

    for j in range(len(data)):
        if j != i:
            training_set.extend(copy(data[j]))
    return (training_set, test_set)

def _return_empty_structure(data):
    if type(data) is CE.Configurations:
        return CE.Configurations()
    elif type(data) is CE.Clusters:
        return CE.Clusters()
    elif type(data) is list:
        return []
