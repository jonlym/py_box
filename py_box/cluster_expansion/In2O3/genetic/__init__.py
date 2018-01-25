from py_box import get_time, get_RMSE, basen_to_base10
from py_box.cluster_expansion import get_correlation_matrix
from sklearn import linear_model
from sklearn.model_selection import train_test_split, cross_val_score, RepeatedKFold
import random
import numpy as np
import pickle
try:
    from mpi4py import MPI
except:
    pass

def save_population(COMM = None, population = None, file_name = None):
	if COMM is None:
		rank = 0
	else:
		rank = COMM.rank

	if rank == 0: 
		with open(file_name, 'w') as f_ptr:
			pickle.dump(population, f_ptr)

def get_bit_string(individual):
	bitstring = ''
	for chromosome in individual:
		if chromosome:
			x = '1'
		else:
			x = '0'
		bitstring = '{}{}'.format(bitstring, x)
	return bitstring

def get_hex_string(individual):
	bitstring = get_bit_string(individual)
	return hex(int(bitstring, 2))

def compare_individuals(individual1, individual2):
    out = np.zeros(len(individual1))
    for i, (x1, x2) in enumerate(zip(individual1, individual2)):
        if x1 == x2:
            out[i] = 1
    return out

def get_similarity(individual1, individual2):
    return np.sum(compare_individuals(individual1, individual2))/len(individual1)

def get_rank(COMM = None):
    if COMM is None:
        return 0
    else:
        return COMM.rank

def get_size(COMM = None):
    if COMM is None:
        return 1
    else:
        return COMM.size

def individual_generator(clusters, n = 0.5):
    for cluster in clusters:
        #Always include the empty cluster
            if cluster.n_nodes == 0:
                yield True
            else:
                rnd_num = random.random()
                yield (rnd_num/cluster.n_nodes > n)

def get_ID(individual):
    return basen_to_base10(individual, 2)    

def get_clusters(individual, clusters_all):
    indices = [i for i, j in enumerate(individual) if j]
    return clusters_all.get_copy(indices = indices)

def split_data(configs_all, E_all, pi_all, kfold):
    indices_train, indices_test, E_train, E_test = train_test_split(range(len(configs_all)), E_all, test_size=1./kfold)
    configs_train = configs_all.get_copy(indices = indices_train)
    configs_test = configs_all.get_copy(indices = indices_test)
    pi_train = pi_all[indices_train]
    pi_test = pi_all[indices_test]
    return (configs_train, configs_test, pi_train, pi_test, E_train, E_test)

def evaluate(individual, clusters_all, configs_all, n_repeats = 1000, kfold = 10):
    clusters = get_clusters(individual, clusters_all)
    E_all = configs_all.get_E_fit(update = True)
    pi_all = get_correlation_matrix(configurations = configs_all, clusters = clusters)

    rmses = np.zeros(shape = (n_repeats, 1))
    rkf = RepeatedKFold(n_splits = kfold, n_repeats = n_repeats)
    regr = linear_model.LinearRegression(fit_intercept = False, n_jobs = -1)
    return (np.abs(np.mean(cross_val_score(estimator = regr, X = pi_all, y = E_all, cv = rkf, scoring = 'neg_mean_squared_error'))),)

    # for i in xrange(n_repeats):
    #     (configs_train, configs_test, pi_train, pi_test, E_train, E_test) = split_data(configs_all = configs_all, E_all = E_all, pi_all = pi_all, kfold = kfold)
    #     regr = linear_model.LinearRegression(fit_intercept = False, n_jobs = -1)
    #     regr.fit(pi_train, E_train)
    #     E_test = regr.predict(pi_test)
    #     rmse = get_RMSE(configs_test.get_E_fit(), E_test)
    #     rmses[i] = rmse
    # return (np.mean(rmses),)

def make_initial_population(COMM = None, toolbox = None, n = None):
    rank = get_rank(COMM)
    if rank == 0:
        print '\t{}  Core {}  Building initial population'.format(get_time(), rank)
        population = toolbox.population(n = n)
    else:
        population = None
    return population


def evaluate_population(COMM = None, toolbox = None, population = None):
    rank = get_rank(COMM)
    size = get_size(COMM)
    if rank == 0:
            jobs_split = np.array_split(range(len(population)), size)
            population_split = []
            for jobs in jobs_split:
                x = []
                for job in jobs:
                    x.append(population[job])
                population_split.append(x)
            print '\t{}  Core {}  Distributing individuals'.format(get_time(), rank)
    else:
        population_split = None
        jobs_split = None

    if COMM is None:
        population_mpi = population
        jobs_mpi = range(len(population))
    else:
        population_mpi = COMM.scatter(population_split, root = 0)
        jobs_mpi = COMM.scatter(jobs_split, root = 0)

    #Evaluate fitness
    fitnesses_mpi = {}
    for i, individual_mpi in zip(jobs_mpi, population_mpi):
        fitnesses_mpi[i] = toolbox.evaluate(individual_mpi)
    print '\t{}  Core {}  Finished evaluating individuals'.format(get_time(), rank)
    if COMM is None:
        fitnesses_list = [fitnesses_mpi]
    else:
        fitnesses_list = MPI.COMM_WORLD.gather(fitnesses_mpi, root = 0)
    if rank == 0:
        print '\t{}  Core {}  Assigning fitness to population.'.format(get_time(), rank)
        for fitnesses_dict in fitnesses_list:
            for i, fitness in fitnesses_dict.iteritems():
                population[i].fitness.values = fitness
    else:
        population = None
    return population

def generate_offspring(COMM = None, toolbox = None, population = None, cxpb = None):
    rank = get_rank(COMM)
    if rank == 0:
        print '\t{}  Core {}  Generating offspring'.format(get_time(), rank)
        print '# Offspring before selection: {}'.format(len(population))
        offspring = toolbox.select(population)
        print '# Offspring after selection: {}'.format(len(offspring))
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

def mutate_offspring(COMM = None, toolbox = None, population = None):
    rank = get_rank(COMM)
    if rank == 0:
        print '\t{}  Core {}  Mutating offspring'.format(get_time(), rank)
        #print '\tApplying mutations'
        for mutant in population:
            toolbox.mutate(mutant)
            del mutant.fitness.values
    else:
        population = None
    return population

def make_next_population(COMM = None, population = None, offspring = None):
    rank = get_rank(COMM)
    if rank == 0:
        print '\t{}  Core {}  Generating new offspring'.format(get_time(), rank)
        population[:] = offspring
    return population

def calculate_statistics(COMM = None, generation = None, population = None):
    rank = get_rank(COMM)
    if rank == 0:
        print '\t{}  Core {}  Calculating statistics'.format(get_time(), rank)
        fitnesses = get_fitnesses(population)
        avg = np.mean(fitnesses)
        sd = np.std(fitnesses)
        min_val = np.min(fitnesses)
        max_val = np.max(fitnesses)
        with open('stats.out', 'a') as f_ptr:
            f_ptr.write('{}\t{}\t{}\t{}\t{}\n'.format(generation, avg, sd, min_val, max_val))

def print_generation_number(COMM = None, generation = None):
    rank = get_rank(COMM)
    if rank == 0:
        print '{}  Core {}  Generation {}'.format(get_time(), rank, generation)

def find_best_individual(COMM = None, population = None):
    rank = get_rank(COMM)
    if rank == 0:
        fitnesses = get_fitnesses(population)
        i = np.where(fitnesses == min(fitnesses))[0][0]
        print '\tIndividual with best fitness:'
        print '\tFitness = {} eV^2'.format(population[i].fitness.values[0])
        print '\tCV RMSE = {} eV'.format(np.sqrt(population[i].fitness.values[0]))
        print '\tHex code: {}'.format(get_hex_string(population[i]))

def get_fitnesses(population):
    return np.array([individual.fitness.values for individual in population])
