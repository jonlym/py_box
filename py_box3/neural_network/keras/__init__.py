import numpy as np
import tensorflow as tf
from keras import backend as K
from keras import initializers, activations, regularizers, constraints
from keras.models import Sequential
from keras.engine.topology import Layer
from keras.layers import Dense, Concatenate

class In2O3Layer(Layer):
	def __init__(self, 
				 activation = None,
				 kernel_initializer = None,
				 bias_initializer = None,
				 kernel_regularizer = None,
				 bias_regularizer = None,
				 activity_regularizer = None,
				 kernel_constraint = None,
				 bias_constraint = None,
				 **kwargs):
		self.activation = activations.get(activation)
		self.kernel_initializer = initializers.get(kernel_initializer)
		self.bias_initializer = initializers.get(bias_initializer)
		self.kernel_regularizer = regularizers.get(kernel_regularizer)
		self.bias_regularizer = regularizers.get(bias_regularizer)
		self.activity_regularizer = regularizers.get(activity_regularizer)
		self.kernel_constraint = constraints.get(kernel_constraint)
		self.bias_constraint = constraints.get(bias_constraint)		
		self.flip_mat = K.constant(value = np.array([
			[0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
			[0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
			[0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.],
			[0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.],
			[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.],
			[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.],
			[1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
			[0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
			[0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
			[0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.],
			[0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
			[0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.]]))
		super(In2O3Layer, self).__init__(**kwargs)

	def build(self, input_shape = None, n_O = 6, n_chain = 2):
		input_shape = n_O * n_chain

		self.kernel = self.add_weight(
			shape=((n_O * n_chain, n_O)), 
			initializer = self.kernel_initializer,
			trainable = True,
			name = 'kernel',
			regularizer = self.kernel_regularizer,
			constraint = self.kernel_constraint)

		self.bias = self.add_weight(
			shape = ((n_O,)),
			initializer = self.bias_initializer,
			trainable = True,
			name = 'bias',
			regularizer = self.bias_regularizer,
			constraint = self.bias_constraint)

		super(In2O3Layer, self).build(input_shape = input_shape)
		self.built = True

	def call(self, inputs):
		output1 = K.bias_add(K.dot(inputs, self.kernel), self.bias)
		output2 = K.bias_add(K.dot(self._flip_input(inputs), self.kernel), self.bias)		
		output = Concatenate(axis = 1)([output1, output2])
		if self.activation is not None:
			output = self.activation(output)
		return output

	def compute_output_shape(self, input_shape):
		return (12, 12) #SHOULD NOT BE HARD CODED

	def _flip_input(self, inputs, n_O = 6):
		return K.dot(inputs, self.flip_mat)


def add_conv_layer(hyperparameters, model, i = 0):
	#Delete the input shape parameter if not the first layer
	if i > 0:
		hyperparameters.pop('input_shape', None)

	model.add(In2O3Layer(**hyperparameters['l1']))
	model.add(Dense(**hyperparameters['l2']))
	return model

def get_hp_cv(hyperparameters, X, y):
	model = build_model(hyperparameters = hyperparameters)
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = hyperparameters['valid_size'])

