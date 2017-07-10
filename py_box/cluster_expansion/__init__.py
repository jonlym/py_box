import numpy as np
from py_box.cluster_expansion.cluster import Cluster
from py_box.cluster_expansion.clusters import Clusters
from py_box.cluster_expansion.configuration import Configuration
from py_box.cluster_expansion.configurations import Configurations
from warnings import warn

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
            interaction_prod = 1
            for index in interaction:
                interaction_prod *= configuration.sigma[index]
            interactions_sum += interaction_prod
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

def get_energies(correlation_mat, Js, configurations = None):
    Es = np.dot(correlation_mat, Js)
    if configurations is not None:
        for E, configuration in zip(Es, configurations):
            configuration.E_CE = E
    return Es

def get_CV(xs_data, xs_fit):
    return np.sqrt(np.mean([(x_data-x_fit)**2 for x_data, x_fit in zip(xs_data, xs_fit)]))
