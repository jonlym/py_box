# coding=utf-8
from py_box3.cluster_expansion.cluster import Cluster
import numpy as np
from copy import copy
from warnings import warn

class Clusters(object):

    def __init__(self, clusters, info = ''):
        self._clusters = list(clusters)
        self.info = info

    def append(self, cluster):
        self._clusters.append(cluster)

    def extend(self, clusters):
        self._clusters.extend(clusters)

    def index(self, name):
        for i, cluster in enumerate(self._clusters):
            if cluster.name == name:
                return i
        else:
            return None

    def remove(self, name = None, index = None):
        if name is not None and index is not None:
            raise Exception('Both name and index cannot be specified.')
        elif name is None and index is None:
            raise Exception('Either name or index must be specified.')
        elif name is not None:
            index = self.index(name)
        del self._clusters[index]

    def get_cluster_energy(self):
        return np.array([cluster.J for cluster in self])

    def get_copy(self, indices = []):
        """Creates a new Clusters object that will contain the indices requested.
        If no indices are specified, then copy the whole object."""
        clusters = []
        #If indices not specified, copy whole object
        if len(indices) == 0:
            indices = list(range(len(self)))

        for i in indices:
            clusters.append(copy(self[i]))
        return Clusters(clusters = clusters)

    def set_Js(self, Js):
        """Assigns Js to the clusters."""
        for J, cluster in zip(Js, self):
            cluster.J = J

    def print_nonsparse(self):
        """Prints the non-sparse interaction terms."""
        for cluster in self:
            if cluster.J != 0.:
                print(('{}\t{}'.format(cluster.name, cluster.J)))

    def print_sparse(self):
        for cluster in self:
            if cluster.J == 0.:
                print(('{}\t{}'.format(cluster.name, cluster.J)))

    def get_n_nodes(self):
        return [cluster.n_nodes for cluster in self]

    def get_domain_matrix(self, concentration):
        """
        Returns the domain matrix described in:
        Mueller, T.; Ceder, G. Phys. Rev. B - Condens. Matter Mater. Phys. 2010, 82 (18), 1–6.
        """
        n_nodes = self.get_n_nodes()
        n_alpha, n_beta = np.meshgrid(n_nodes, n_nodes)
        return (2*concentration - 1)**(n_alpha + n_beta)

    def __len__(self):
        return len(self._clusters)

    def __setitem__(self, index, cluster):
        self._clusters[index] = cluster

    def __getitem__(self, index):
        return self._clusters[index]

    @classmethod
    def from_excel(cls, file_name = 'clusters.xlsx', skiprows = [], info = ''):
        import pandas as pd

        if info == '':
            info = 'Data generated using class method, from_excel(), using file: {}. Rows skipped: {}'.format(file_name, skiprows)
        cluster_data = pd.read_excel(file_name, skiprows = skiprows, header=0)
        clusters = []
        for i, cluster in cluster_data.iterrows():
            for col in cluster_data.columns:
                #Process blank records
                if type(cluster[col]) is float and np.isnan(cluster[col]):
                    cluster.set_value(col, None)
            #Convert string into list of list
            out_interactions = []
            if cluster['interactions'] is not None:
                if '|' in cluster['interactions']:
                    interactions = cluster['interactions'].split('|')
                else:
                    interactions = [cluster['interactions']]

                for interaction in interactions:
                    if type(interaction) is int:
                        out_interactions.append([interaction])
                    elif type(interaction) is float:
                        out_interactions.append([int(interaction)])
                    elif type(interaction) is str:
                        interaction_indices = []
                        for interaction_index in interaction.split(','):
                            interaction_indices.append(int(interaction_index))
                        out_interactions.append(interaction_indices)
                    else:
                        warn('Invalid data type for record {}: {}'.format(cluster['name'], type(interaction)))

            clusters.append(Cluster(name = cluster['name'],
                                    m = cluster['m'],
                                    N = cluster['N'],
                                    n_cell = cluster['n_cell'],
                                    interactions = out_interactions,
                                    info = cluster['info'],
                                    J = cluster['J'],
                                    n_nodes = cluster['n_nodes']))
        return cls(clusters = clusters, info = info)
