class Cluster(object):
    """Describes the cluster.
        Attributes:
            name - string
                Name of cluster. Required as a search key
            m - int
                Holds the multiplicity (or degeneracy) of the cluster in the considered unit cell.
            N - int
                Number of vertices
            n_cell - int
                Size of unit cell
            interactions - List of lists of int
                Indices of the configuration vector to test
                Each list contains the product that will be taken.
            J - float
                Effective cluster interaction (eV)
            n_nodes - int
                Number of nodes in the cluster
            info - string
                Description of the cluster

    """
    type_dict = {'name': str,
                 'm': int,
                 'N': int,
                 'n_cell': int,
                 'n_nodes': int,
                 'interactions': list,
                 'J': float,
                 'info': str}

    def __init__(self, name = None, m = None, N = None, n_cell = None, interactions = None, J = None, n_nodes = None, info = None):
        if name is not None:
            self.name = str(name)
        else:
            self.name = None
        self.m = m
        self.N = N
        self.n_cell = n_cell
        self.interactions = interactions
        self.J = J
        self.n_nodes = n_nodes
        if info is not None:
            self.info = str(info)
        else:
            self.info = None

