import numpy as np
import networkx as nx

z_O = 18.196
fcc_sites = {0:  np.array([1.403, 0.810, z_O]),
             1:  np.array([2.805, 3.239, z_O]),
             2:  np.array([4.208, 5.669, z_O]),
             3:  np.array([5.611, 8.098, z_O]),
             4:  np.array([4.208, 0.810, z_O]),
             5:  np.array([5.611, 3.329, z_O]),
             6:  np.array([7.013, 5.669, z_O]),
             7:  np.array([8.416, 8.098, z_O]),
             8:  np.array([7.013, 0.810, z_O]),
             9:  np.array([8.416, 3.239, z_O]),
             10: np.array([9.819, 5.669, z_O]),
             11: np.array([11.221, 8.098, z_O]),
             12: np.array([9.819, 0.810, z_O]),
             13: np.array([11.221, 3.239, z_O]),
             14: np.array([12.624, 5.669, z_O]),
             15: np.array([14.027, 8.098, z_O]),
             }

def get_cell(size = 4., a = 1.):
    size = float(size)
    a = float(a)

    x_offset = a/2.
    y_offset = np.sqrt(3.)/2. * a
    return np.array([[size*a, 0.], [size*x_offset, size*y_offset]])

def get_sigma(n, size = 4):
    return [int(i) for i in np.binary_repr(n, size**2)]

def get_distance(pos_i, pos_j):
    return np.sqrt(np.sum((i-j)**2 for i, j in zip(pos_i, pos_j)))

def get_distances(pos_i_dict, pos_j_dict):
    distances = {}
    for key_i, val_i in pos_i_dict.iteritems():
        for key_j, val_j in pos_j_dict.iteritems():
            distances[(key_i, key_j)] = get_distance(val_i, val_j)
    return distances

def get_min_distance(pos_i_dict, pos_j_dict):
    return min([val for key, val in get_distances(pos_i_dict, pos_j_dict).iteritems()])
    # distances = []
    # for key_i, val_i in pos_i_dict.iteritems():
    #     for key_j, val_j in pos_j_dict.iteritems():
    #         distances.append(get_distance(val_i, val_j))
    # return min(distances)

def get_big_graph(size = 4, a = 1, eps = None, periodic = True):
    if eps is None:
        eps = a/100.
    big_graph = nx.Graph()
    pos_draw = {}


    #Create cell dimensions
    x_offset = a/2.
    y_offset = np.sqrt(3.)/2. * a

    cell = get_cell(size = size, a = a)

    #Create nodes
    k = 0
    for i in xrange(size):
        for j in xrange(size):
            #Position within cell
            orig_pos = (x_offset*j + a*i, y_offset*j)
            #Find mirror copies of images outside cell
            positions = {}
            if periodic:
                for l in xrange(-1, 2):
                    for m in xrange(-1, 2):
                        direction = np.array([[l], [m]])
                        positions[(l, m)] = np.sum(cell * direction, axis = 0) + orig_pos
            else:
                positions[(0, 0)] = orig_pos
            big_graph.add_node(k, positions = positions)
            k = k + 1
    #Create edges
    for i, data_i in big_graph.nodes_iter(data = True):
        for j, data_j in big_graph.nodes_iter(data = True):
            if i != j:
                min_distance = get_min_distance(data_i['positions'], data_j['positions'])
                if min_distance - a < eps:
                    big_graph.add_edge(i, j, min_distance = min_distance)
    return big_graph

def get_small_graph(config, size = 4, a = 1, eps = None, periodic = True, dict_edges = False):
    big_graph = get_big_graph(size = size, a = a, eps = eps, periodic = periodic)
    #Remake graph
    graph_config = get_sigma(n = config, size = size)
    graph = nx.Graph()
    positions = {}
    pos_draw = {}
    for j, k in enumerate(graph_config):
        if k == 1:
            positions = nx.get_node_attributes(big_graph, 'positions')[j]
            graph.add_node(j, positions = positions)
            pos_draw[j] = positions[(0, 0)]
    #Create edges in the graph
    for (node_i, data_i) in graph.nodes_iter(data = True):
        for (node_j, data_j) in graph.nodes_iter(data = True):
            if node_i != node_j:
                if periodic:
                    graph.add_edge(node_i, node_j, distance = get_min_distance(data_i['positions'], data_j['positions']))
                else:
                    graph.add_edge(node_i, node_j, distance = get_distance(data_i['positions'][(0, 0)], data_j['positions'][(0, 0)]))

                if dict_edges:
                    graph.add_edge(node_i, node_j, distances = get_distances(data_i['positions'], data_j['positions']))
    return graph

def draw_graph_node(config, size = 4, a = 1, eps = None):
    import matplotlib.pyplot as plt

    big_graph = get_big_graph(size = size, a = a, eps = eps)
    graph = get_small_graph(config, size = size, a = a, eps = eps)
    graph_config = [int(j) for j in np.binary_repr(config, size**2)]

    pos_draw = {}
    for j, k in enumerate(graph_config):
        if k == 1:
            pos_draw[j] = nx.get_node_attributes(big_graph, 'positions')[j][(0, 0)]
    plt.figure(config)
    nx.draw(graph, pos = pos_draw, node_size = 700)
    nx.draw_networkx_labels(graph, pos = pos_draw)
    plt.draw()
