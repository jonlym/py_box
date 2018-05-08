import numpy as np

def add_periodic_padding(X, pad_size):
	"""
	Pads the matrix, X, with a periodic boundary condition. 
	The pad size should be passed as an iterable object that should have the same shape as
	the shape of X. i.e. If X is 2 dimensional, pad_size should contain two elements.
	"""
	pad_size = np.array(pad_size)
	n_duplicates = tuple([int(x) for x in np.ceil(pad_size/np.array(X.shape))*2 + 1])
	X_out = np.tile(X, n_duplicates)
	n_dlt = [int(x) for x in (np.array(X.shape) - np.mod(pad_size, np.array(X.shape)))]
	X_out = X_out[:-n_dlt[0], :]
	X_out = X_out[:, :-n_dlt[1]]
	X_out = X_out[n_dlt[0]:, :]
	X_out = X_out[:, n_dlt[1]:]
	return X_out