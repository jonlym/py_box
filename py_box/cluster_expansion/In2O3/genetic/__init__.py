from py_box import get_time, get_RMSE
from py_box.cluster_expansion import get_correlation_matrix
import random
import numpy as np
from sklearn import linear_model
from sklearn.model_selection import train_test_split
try:
	from mpi4py import MPI
except:
	pass

def individual_generator(clusters):
	while True:
	    for cluster in clusters:
	    	#Always include the empty cluster
	    	if cluster.n_nodes == 0:
	    		yield True
	    	else:
		        yield random.random()/cluster.n_nodes > 0.5

def get_clusters(individual, clusters_all):
	indices = [i for i, j in enumerate(individual) if j]
	return clusters_all.get_copy(indices = indices)

def split_data(configs_all, E_all, pi_all, kfold):
	indices_train, indices_test, E_train, E_test = train_test_split(range(len(configs_all)), E_all, test_size=1./kfold)
	configs_test = configs_all.get_copy(indices = indices_test)
	pi_train = pi_all[indices_train]
	pi_test = pi_all[indices_test]
	return (configs_test, pi_train, pi_test, E_train, E_test)

def evaluate(individual, clusters_all, configs_all, n_repeats = 1000, kfold = 10):
	clusters = get_clusters(individual, clusters_all)
        print 'Size of individual: {}'.format(np.sum([1 for x in individual if x]))
        print 'Size of clusters: {}'.format(len(clusters))
        for cluster in clusters:
		print cluster.name
	E_all = configs_all.get_E_fit(update = True)
	pi_all = get_correlation_matrix(configurations = configs_all, clusters = clusters)

	rmses = np.zeros(shape = (n_repeats, 1))
	for i in xrange(n_repeats):
		(configs_test, pi_train, pi_test, E_train, E_test) = split_data(configs_all = configs_all, E_all = E_all, pi_all = pi_all, kfold = kfold)
		regr = linear_model.LinearRegression(fit_intercept = False, n_jobs = -1)
		regr.fit(pi_train, E_train)
		E_test = regr.predict(pi_test)
		rmse = get_RMSE(configs_test.get_E_fit(), E_test)
		rmses[i] = rmse
	return (np.mean(rmses),)

def make_initial_population(COMM, toolbox, n):
	if COMM.rank == 0:
		print '\t{}  Core {}  Building initial population'.format(get_time(), COMM.rank)
		population = toolbox.population(n = n)
	else:
		population = None
	return population


def evaluate_population(COMM, toolbox, population):
	if COMM.rank == 0:
    		jobs_split = np.array_split(range(len(population)), COMM.size)
    		population_split = []
    		for jobs in jobs_split:
    			x = []
    			for job in jobs:
        			x.append(population[job])
    			population_split.append(x)
		print '\t{}  Core {}  Distributing individuals'.format(get_time(), COMM.rank)
	else:
		population_split = None
                jobs_split = None
	population_mpi = COMM.scatter(population_split, root = 0)
	jobs_mpi = COMM.scatter(jobs_split, root = 0)

	#Evaluate fitness
	fitnesses_mpi = {}
	for i, individual_mpi in zip(jobs_mpi, population_mpi):
		fitnesses_mpi[i] = toolbox.evaluate(individual_mpi)
        print '\t{}  Core {}  Finished evaluating individuals'.format(get_time(), COMM.rank)
        fitnesses_list = MPI.COMM_WORLD.gather(fitnesses_mpi, root = 0)
	if COMM.rank == 0:
		print '\t{}  Core {}  Assigning fitness to population.'.format(get_time(), COMM.rank)
        	for fitnesses_dict in fitnesses_list:
	        	for i, fitness in fitnesses_dict.iteritems():
		        	population[i].fitness.values = fitness
	else:
		population = None
	return population

def generate_offspring(COMM, toolbox, population, cxpb):
	if COMM.rank == 0:
        	print '\t{}  Core {}  Generating offspring'.format(get_time(), COMM.rank)
		offspring = toolbox.select(population)
		offspring = [toolbox.clone(individual) for individual in offspring]

		#Apply crossover and mutation on the offspring
		#print '\tMaking children'
		for child1, child2 in zip(offspring[::2], offspring[1::2]):
			if random.random() < cxpb:
				toolbox.mate(child1, child2)
				del child1.fitness.values
				del child2.fitness.values
	else:
		offspring = None
	return offspring

def mutate_offspring(COMM, toolbox, population):
	if COMM.rank == 0:
    		print '\t{}  Core {}  Mutating offspring'.format(get_time(), COMM.rank)
		#print '\tApplying mutations'
		for mutant in population:
			toolbox.mutate(mutant)
			del mutant.fitness.values
	else:
		population = None
	return population

def make_next_population(COMM, population, offspring):
	if COMM.rank == 0:
		print '\t{}  Core {}  Generating new offspring'.format(get_time(), COMM.rank)
		population[:] = offspring
	return population

def calculate_statistics(COMM, generation, population):
	if COMM.rank == 0:
		print '\t{}  Core {}  Calculating statistics'.format(get_time(), COMM.rank)
		fitnesses = get_fitnesses(population)
		avg = np.mean(fitnesses)
		sd = np.std(fitnesses)
		min_val = np.min(fitnesses)
		max_val = np.max(fitnesses)
		with open('stats.out', 'a') as f_ptr:
			f_ptr.write('{}\t{}\t{}\t{}\t{}\n'.format(generation, avg, sd, min_val, max_val))

def print_generation_number(COMM, generation):
	if COMM.rank == 0:
		print '{}  Core {}  Generation {}'.format(get_time(), COMM.rank, generation)

def find_best_individual(COMM, population):
	if COMM.rank == 0:
		fitnesses = get_fitnesses(population)
        	i = np.where(fitnesses == max(fitnesses))[0][0]
        	print 'Individual with best fitness:'
        	print 'Fitness = {} eV'.format(population[i].fitness.values[0])
        	print population[i]

def get_fitnesses(population):
	return np.array([individual.fitness.values for individual in population])
