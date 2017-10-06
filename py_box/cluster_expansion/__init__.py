import numpy as np
from py_box import get_null, get_time
from py_box.cluster_expansion.cluster import Cluster
from py_box.cluster_expansion.clusters import Clusters
from py_box.cluster_expansion.configuration import Configuration
from py_box.cluster_expansion.configurations import Configurations
from py_box.ase.In2O3 import run_In2O3_configuration
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import normalize
from warnings import warn

def get_configuration_difference(train_mat, new_mat):
    null_train_mat = get_null(train_mat)

    #Normalize matrices
    norm_new_mat = normalize(new_mat, axis = 1)
    norm_null_train_mat = normalize(null_train_mat, axis = 0)

    max_dot_prod = np.amax(np.absolute(np.dot(norm_new_mat, norm_null_train_mat)), axis = 1)
    return max_dot_prod

def get_difference(big_obj, small_obj):
    """
    Gets the difference between a larger object and a smaller object. Uses the name to determine if the records are the
    same. Should be compatible with any objects with the attribute name.
    :param big_obj:
    :param small_obj:
    :return:
    """
    diff_obj = big_obj.get_copy()
    for small_item in small_obj:
        diff_obj.remove(name = small_item.name)
    return diff_obj

def get_correlation_matrix(configurations, clusters):
    """Calculates the correlation matrix"""
    correlation_matrix = np.zeros((len(configurations), len(clusters)))
    for i, configuration in enumerate(configurations):
        for j, cluster in enumerate(clusters):
            correlation_matrix[i, j] = _get_correlation(configuration, cluster)
            if correlation_matrix[i, j] > 1:
                warn('Correlation matrix element ({},{}) greater than 1.'.format(i, j))
            elif correlation_matrix[i, j] < -1:
                warn('Correlation matrix element ({},{}) less than -1.'.format(i, j))
    return correlation_matrix

def _get_correlation(configuration, cluster):
    """Calculates an element of the correlation matrix with a specific configuration and cluster."""
    if len(cluster.interactions) == 0:
        return 1 #Empty cluster has a constant value
    else:
        interactions_sum = 0
        for interaction in cluster.interactions:
            interactions_sum += np.prod([configuration.sigma[i] for i in interaction])
        return interactions_sum/(float(cluster.N) * float(cluster.m))

def get_effective_interactions(correlation_mat, energies, clusters = None):
    # Js = np.linalg.lstsq(correlation_mat, energies)
    #Using the formula: b = (X' X)^-1 X' Y from https://onlinecourses.science.psu.edu/stat501/node/382
    pi = correlation_mat
    pi_T = np.transpose(correlation_mat)
    Js = np.dot(np.dot(np.linalg.inv(np.dot(pi_T, pi)), pi_T), energies)
    if clusters is not None:
        for J, cluster in zip(Js, clusters):
            cluster.J = J
    return Js

def get_effective_interactions_lasso(correlation_mat, energies, alpha = [0.1], clusters = None):
    clf = LassoCV(alphas = alpha, fit_intercept=False, copy_X=True)
    clf.fit(correlation_mat, energies)
    Js = clf.coef_
    if clusters is not None:
        for J, cluster in zip(Js, clusters):
            cluster.J = J
    return (Js, clf.alpha_)

def get_energies(correlation_mat, Js, configurations = None, intercept = 0.):
    Es = np.dot(correlation_mat, Js) + intercept
    if configurations is not None:
        for E, configuration in zip(Es, configurations):
            configuration.E_CE = E
    return Es

def get_RMSE(xs_data, xs_fit):
    return np.sqrt(np.mean([(x_data-x_fit)**2 for x_data, x_fit in zip(xs_data, xs_fit)]))

def count_nonsparse(Js, eps = 1e-5):
    return np.count_nonzero([abs(J) > eps for J in Js])

def get_best_structure(energies, differences, cv, cv_limit, n = 0):
    """
    Returns the index of the 'best' structure to add to the optimization based on the optimization of energies and differences
    :param energies:
    :param differences:
    :return:
    """
    if cv < cv_limit:
        return np.where(energies == sorted(energies, reverse = False)[n])[0][0]
    else:
        return np.where(differences == sorted(differences, reverse = True)[n])[0][0]

def run_cluster_expansion(train_path, clusters_path, configs_all_path, log_path, submit_job = True, n_new = 1):
    #Read cluster data
    print 'Reading cluster data'
    clusters = Clusters.from_excel(clusters_path)

    #Read training structures
    print 'Reading configuration training data'
    configs_train = Configurations.from_vasp(train_path)
    configs_train.set_E_fit()
    #Read all training structures
    print 'Reading all configuration data'
    configs_all = Configurations.from_excel(configs_all_path)
    print 'Finding difference'
    configs_new = get_difference(configs_all, configs_train)

    #Generating correlation matrices
    print 'Generating correlation matrix for training structures'
    pi_train = get_correlation_matrix(configurations = configs_train, clusters = clusters)
    print 'Generating correlation matrix for new structures'
    pi_new = get_correlation_matrix(configurations = configs_new, clusters = clusters)

    #Find structures that would result in better CV score
    print 'Calculating similarity of new configurations to training configurations'
    configs_difference = get_configuration_difference(pi_train, pi_new)

    #Run Cluster Expansion Model
    print 'Running Lasso with Leave-One-Out Cross Validation'
    clf = LassoCV(copy_X=True, cv = len(configs_train), fit_intercept = True)
    print configs_train.get_E_fit()
    clf.fit(pi_train, configs_train.get_E_fit())

    #Print Model Data
    Js = clf.coef_
    intercept = clf.intercept_

    #Calculate energies
    print 'Calculating energies using Cluster Expansion'
    CE_E_new = get_energies(correlation_mat = pi_new, Js = Js, intercept = intercept)

    #Start DFT calculation for structure
    j = 0
    new_structures = []
    new_indices = []
    for n in xrange(len(configs_new)):
        new_index = get_best_structure(CE_E_new, configs_difference, cv = np.average(clf.mse_path_[-1]), cv_limit = 0.0025, n = n)
        print 'Attempting to submit {}'.format(configs_new[new_index].name)
        successful_submit = run_In2O3_configuration(configs_new[new_index], rel_path = train_path, submit_job = submit_job)
        if successful_submit:
            j = j + 1
            new_structures.append(configs_new[new_index].name)
            new_indices.append(new_index)
        else:
            print 'Failed to submit {}'.format(configs_new[new_index].name)

        if j >= n_new:
            break
    else:
        print 'Could not find structure to submit.'
        new_structures = ['Nan']

    print 'Updating log file, {}'.format(log_path)
    with open(log_path, 'a') as log_ptr:
        log_ptr.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(get_time(), clf.alpha_, np.average(clf.mse_path_[-1]), count_nonsparse(Js = Js), '\t'.join([new_structure for new_structure in new_structures]), '\t'.join([str(configs_difference[new_index]) for new_index in new_indices]), '\t'.join([str(CE_E_new[new_index]) for new_index in new_indices])))
