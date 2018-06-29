import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
import math
from py_box3.ase import get_distance
from matplotlib import pyplot as plt

def get_111_graph(r = 1, a = 1.):
	"""
	Generates a (111) graph with a radius, r, and bond lengths, a.

	Parameters
	----------
		r - int
			Radius of the graph
		a - float
			Bond length
	Returns
	-------
		111_graph - networkx.Graph object
			A graph with the connectivity similar to a (111) lattice.
	"""

	G = nx.Graph()
	previous_positions = []
	new_layers_indices = [0, 1]
	for i in range(r):
		new_layers_indices.append(new_layers_indices[-1] + 6 * (len(new_layers_indices)-1))
	n = new_layers_indices[-1]
	for i in range(n):
		#First node centered at origin
		if i in new_layers_indices:
			new_layer = True
		else:
			new_layer = False

		if i == 0:
			layer = 0
			position = np.array([0., 0.])
		#Second node displayed from origin by a
		elif i == 1:
			layer = 1
			position = np.array([a, 0.])
		#Positions of subsequent nodes inferred from previous nodes
		else:
			pos1 = G.nodes[i - 1]['position']
			#Creating a new layer
			if new_layer:
				layer += 1
				pos2 = G.nodes[new_layers_indices[layer-1]]['position']
			else:
				#Find node connected to i-1
				if (i-1 in new_layers_indices) and (i != 0):
					j = min([x for x in G.neighbors(i-1) if G.nodes[x]['layer'] == layer-1])
				else:
					j = max([x for x in G.neighbors(i-1) if G.nodes[x]['layer'] == layer-1])
				pos2 = G.nodes[j]['position']
			new_positions = get_perpendicular_bisectors(position1 = pos1, position2 = pos2, a = a)
			position = get_unique_position(new_positions = new_positions, previous_positions = previous_positions)

		node_attr = {
			'layer': layer,
			'position': position
			}
		G.add_node(i, **node_attr)
		G = connect_new_edges(new_node = i, graph = G, new_layer = new_layer, a = a)
		previous_positions.append(position)
	return G

def get_perpendicular_bisectors(position1, position2, a):
	"""
	Calculates the location of the perpendicular bisectors from two points.
	The output positions are 'a' away from the input positions.

	Parameters
	----------
		position1 - (2,) ndarray
			x and y coordinate for point 1
		position2 - (2,) ndarray
			x and y coordinates for point 2
		a - float
			Distance away from two points
	Returns
	-------
		positions_out - (2, 2) ndarray
			Coordinates of two perpendicular bisectors
	"""
	mid_point = np.mean(np.vstack([position1, position2]), axis = 0)
	v_old = position1 - position2
	if math.isclose(v_old[1], 0.):
		v_new = np.array([0., np.sqrt(3)/2 * a])
	else:
		v_new = np.array([0., 0.])
		v_new[0] = np.sqrt(0.75 * a**2 * 1/(1 + (v_old[0]/v_old[1])**2))
		v_new[1] = -v_old[0] * v_new[0] / v_old[1]
	positions_out = np.array([mid_point + v_new, mid_point - v_new])
	return positions_out

def get_unique_position(new_positions, previous_positions):
	"""
	Returns the new location given a list of previous positions.
	
	Parameters
	----------
		new_positions - (X, Y) ndarray
			Possible options for the new position
		previous_positions - (Z, Y) ndarray
			Previous positions found
	Returns
	-------
		unique_positions - (Y,) ndarray
			The unique position not found in previous_position 
	"""
	for new_position in new_positions:
		if not np.any([np.allclose(new_position, x) for x in previous_positions]):
			return new_position

def connect_new_edges(new_node, graph, new_layer, a):
	"""
	Connects the new edges to the graph.

	Parameters
	----------
		new_node - int
			Name of the new node
		graph - networkx.Graph object
			Graph where connections will be made
	Returns
	-------
		graph - networkx.Graph object
			Updated graph object
	"""
	position = graph.nodes[new_node]['position']
	layer = graph.nodes[new_node]['layer']

	#New node is connected to last index. New layers will be connected by layer rule
	if not new_layer:
		graph.add_edge(new_node, new_node-1)
	for node in graph.nodes(data = True):
		#Nodes are not connected to themselves
		if node[0] == new_node:
			continue
		#New node is only connected to close-by nodes
		distance = get_distance(position, node[1]['position'])
		if not math.isclose(a, distance):
			continue
		#If criteria met, add edge
		graph.add_edge(new_node, node[0])
	return graph

def get_unoccupied_nodes(subgraph, graph):
	"""
	Finds the nodes in large_graph that are not present in graph

	Parameters
	----------
		subgraph - nx.Graph()
			The smaller graph that is missing nodes from large_graph
		graph - nx.Graph()
			Larger graph that contains all the relevant nodes
	Returns
	-------
		unoccupied_nodes - (N,) tuple
			Nodes present in graph but not present in subgraph
	"""

	subgraph_nodes = set(subgraph.nodes())
	graph_nodes = set(graph.nodes())

	return graph_nodes.difference(subgraph_nodes)

def get_bit_string(graph, n):
	"""
	Returns the bitstring representation of the graph.

	Parameters
		graph - networkx.Graph object
			The graph to encode. The node names should be ascending integers
		n - int
			Length of desired bitstring
	Returns
		str
			The bitstring representation of the graph
	"""

	bitstring_rev = ''
	for i in range(n):
		if graph.has_node(i):
			new_character = 1
		else:
			new_character = 0
		bitstring_rev = '{}{}'.format(bitstring_rev, new_character)
	return bitstring_rev[::-1]

def get_isomorphic_graphs(large_graph, small_graph):
	"""
	Returns the isomorphic graphs
	"""
	GM = iso.GraphMatcher(large_graph, small_graph, edge_match=iso.numerical_edge_match(attr = 'distance', default = 1.))
	isomorphic_dicts = [x for x in GM.subgraph_isomorphisms_iter()]
	unique_isomorphic_sets = get_unique_isomorphic_graphs(isomorphic_dicts)
	return unique_isomorphic_sets

def get_unique_isomorphic_graphs(isomorphic_dicts):
	"""
	Eliminates redundant isomorphic graphs.

	Parameters
	----------
		isomorphic_dicts - dict
			Dictionary generated using networkx.algorithms.isomorphism.GraphMatcher
	Returns
	-------
		unique_isomorphic_set - set
			Set of sets with nodes
	"""
	return set([frozenset(isomorphic_dict.keys()) for isomorphic_dict in isomorphic_dicts])