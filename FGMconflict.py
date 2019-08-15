#!/usr/bin/env python

### This program simulates two populations evolving under Fisher's geometric model with conflict and a control without conflict ###
### python3 FGMconflict.py -help for input options ###
### Written by Trey J Scott 2018 ###
### python --version ###
### Python 3.5.2 :: Anaconda 4.2.0 (x86_64) ###

# Import programs
import random
import numpy as np
from scipy.spatial import distance as dist
from scipy.stats import norm
import scipy.stats as stats
import argparse
import scipy.special as spc

### FUNCTIONS ###

# Function to generate random mutations with a specified average size
def generate_random_vector(average_mutation):
	sd_1d = average_mutation*((np.pi)**(1/2))/(2**(1/2))
	uni = 2*average_mutation
	expo = average_mutation
	if distribution == 'uniform':
		radial = np.random.uniform(0,uni)
	if distribution == 'exponential':
		radial = np.random.exponential(expo)
	if distribution == 'normal':
		radial = abs(np.random.normal(0, sd_1d))
		vector = np.array(radial * (-1)**random.randint(1,2))
	return radial, vector
# Generates optima for both parties
def generate_optima(d1, d2):
	optimum1 = np.array([(-(1/d1)*np.log(r))**(1/Q)])
	optimum2 = np.array([-(-(1/d2)*np.log(r))**(1/Q)])
	return optimum1, optimum2

# Gaussian fitness function
def fitness_function(distance, d):
	return np.exp(-(d*(distance**Q)))

# Calculates probability of fixation for new mutations
def calculate_u(new_distance, old_distance, N = 'infinite', denominator = 2):
	fitness_new = fitness_function(new_distance, denominator)
	fitness_old = fitness_function(old_distance, denominator)
	s_coefficient = (fitness_new/fitness_old) - 1
	if N == 'infinite':
		probability_of_fixation = (1 - np.exp(-2*s_coefficient))
	elif N > 0:
		probability_of_fixation = ((1 - np.exp(-2*s_coefficient))/(1 - np.exp(-4*s_coefficient*N)))
	return probability_of_fixation, s_coefficient

# Function that simulates standard adaptation with Fisher's geometric model
def standard_adaptation(position, optimum, samp):
	counter = 0
	distance_to_optimum = dist.euclidean(position, optimum)
	while counter < m:
		mutation_size, vector = generate_random_vector(average_mutation1)
		future_position = position + vector
		new_dist_to_optimum = dist.euclidean(future_position, optimum)
		u, s = calculate_u(new_dist_to_optimum, distance_to_optimum, N_1, d1)
		if random.random() <= u:
			mutation_fitness = vector
			position = future_position
			distance_to_optimum = new_dist_to_optimum
			if counter >= burn_in:
				output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s) + ',' + str(mutation_size) + ',' + str(fitness_function(distance_to_optimum,d1)) + ',No Conflict,No Conflict,Fixed\n')
		else:
			if counter >= burn_in:
				output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s) + ',' + str(mutation_size) + ',' + str(fitness_function(distance_to_optimum,d1)) + ',No Conflict,No Conflict,Unfixed\n')
		counter += 1

# Function that simulates conflict over a joint phenotype with Fisher's geometric model
def conflict_model(position, optimum1, optimum2, samp):
	position = position
	counter = 0
	distance_to_optimum1 = dist.euclidean(position, optimum1)
	distance_to_optimum2 = dist.euclidean(position, optimum2)
	while counter < m:
		# Test which party mutates first
		if random.randint(1,2) == 1:
			# party 1 will mutate j times for every 1 mutation of party 2
			for j in range(rate):
				mutation_size1, random_vector = generate_random_vector(average_mutation1)
				future_position = position + random_vector
				new_dist_to_optimum1 = dist.euclidean(future_position, optimum1)
				u1, s1 = calculate_u(new_dist_to_optimum1, distance_to_optimum1, N_1, d1)
				if random.random() <= u1:
					position = future_position
					distance_to_optimum1 = new_dist_to_optimum1
					distance_to_optimum2 = dist.euclidean(position, optimum2)
					mut_out.write(str(random_vector) + ',Population 1\n')
					if counter >= burn_in:
						output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s1) + ',' + str(mutation_size1) + ',' + str(fitness_function(distance_to_optimum1,d1)) + ',Conflict,Party 1,Fixed\n')
				else:
					if counter >= burn_in:
						output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s1) + ',' + str(mutation_size1) + ',' + str(fitness_function(distance_to_optimum1,d1)) + ',Conflict,Party 1,Unfixed\n')
					mut_out.write('0,Population 1\n')
			### Party 2 mutation
			mutation_size2, random_vector = generate_random_vector(average_mutation2)
			future_position = position + random_vector
			new_dist_to_optimum2 = dist.euclidean(future_position, optimum2)
			u2, s2 = calculate_u(new_dist_to_optimum2, distance_to_optimum2, N_2, d2)
			if random.random() <= u2:
				position = future_position
				distance_to_optimum1 = dist.euclidean(position, optimum1)
				distance_to_optimum2 = new_dist_to_optimum2
				mut_out.write(str(random_vector) + ',Population 2\n')
				if counter >= burn_in:
					output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s2) + ',' + str(mutation_size2) + ',' + str(fitness_function(distance_to_optimum2,d2)) + ',Conflict,Party 2,Fixed\n')
			else:
				if counter >= burn_in:
					output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s2) + ',' + str(mutation_size2) + ',' + str(fitness_function(distance_to_optimum2,d2)) + ',Conflict,Party 2,Unfixed\n')
				mut_out.write('0,Population 2\n')
			counter += 1
		else:
			### Party 2 mutates
			mutation_size2, random_vector = generate_random_vector(average_mutation2)	
			future_position = position + random_vector
			new_dist_to_optimum2 = dist.euclidean(future_position, optimum2)
			u2, s2 = calculate_u(new_dist_to_optimum2, distance_to_optimum2, N_2, d2)
			if random.random() <= u2:
				position = future_position
				distance_to_optimum1 = dist.euclidean(position, optimum1)
				distance_to_optimum2 = new_dist_to_optimum2
				mut_out.write(str(random_vector) + ',Population 2\n')
				if counter >= burn_in:
					output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s2) + ',' + str(mutation_size2) + ',' + str(fitness_function(distance_to_optimum2,d2)) + ',Conflict,Party 2,Fixed\n')
			else:
				if counter >= burn_in:
					output.write(str(counter) + ',' + str(samp)+ ',' + str(position[0]) + ',' + str(s2) + ',' + str(mutation_size2) + ',' + str(fitness_function(distance_to_optimum2,d2)) + ',Conflict,Party 2,Unfixed\n')
				mut_out.write('0,Population 2\n')
			### Party 1 mutates
			for j in range(rate):
				mutation_size1, random_vector = generate_random_vector(average_mutation1) 
				future_position = position + random_vector
				new_dist_to_optimum1 = dist.euclidean(future_position, optimum1)
				u1, s1 = calculate_u(new_dist_to_optimum1, distance_to_optimum1, N_1, d1)
				if random.random() <= u1:
					position = future_position
					distance_to_optimum1 = new_dist_to_optimum1
					distance_to_optimum2 = dist.euclidean(position, optimum2)
					mut_out.write(str(random_vector) + ',Population 1\n')
					if counter >= burn_in:
						output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s1) + ',' + str(mutation_size1) + ',' + str(fitness_function(distance_to_optimum1,d1)) + ',Conflict,Party 1,Fixed\n')
				else:
					if counter >= burn_in:
						output.write(str(counter) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s1) + ',' + str(mutation_size1)+ ',' + str(fitness_function(distance_to_optimum1,d1)) + ',Conflict,Party 1,Unfixed\n')
					mut_out.write('0,Population 1\n')
			counter += 1

# Runs multiple simulations
def run_simulation(position, num_samples):
	optimum1, optimum2 = generate_optima(d1,d2)
	for sample in range(num_samples):
		standard_adaptation(position, optimum1, sample)
		conflict_model(position, optimum1, optimum2, sample)
	output.close()
	mut_out.close()

### SET ARGUMENTS
ap = argparse.ArgumentParser()
ap.add_argument('-x', '--samples', help = 'number of replicate simulations (default is 500)', type = int)
ap.add_argument('-p', '--population_size1', help = 'population size for first party (default is infinite)', type = int)
ap.add_argument('-pp', '--population_size2', help = 'population size for second party (defualt is infinite)', type = int)
ap.add_argument('-m', '--mutations', help = 'mutation distribution for mutation vectors (normal-default, uniform, or exponential)')
ap.add_argument('-q', '--Q', help = 'changes Q (epistasis) parameter in fitness function', type = float)
ap.add_argument('-z', '--attempts', help = 'number of iterations per single simulation (default is 5000)', type = int)
ap.add_argument('-c', '--init_fit', help = 'defines initial fitness (default is 0.2)', type = float)
ap.add_argument('-r', '--rate', help = 'allows first party to mutate r times for every one mutation of second party', type = int)
ap.add_argument('-b', '--burn_in', help = 'defines burn in period for equilibrium', type = int)
ap.add_argument('-a', '--ave_mut', help = 'scales average mutation sizes for both parties (default is 0.2)', type = float)
ap.add_argument('-aa', '--rel_ave_mut', help = 'scales average mutation sizes for Party 1 relative to Party 2', type = float)
ap.add_argument('-d', '--selection', help = 'Adjust strength of selection for both parties', type = float)
ap.add_argument('-dd', '--rel_selection', help = 'Adjust strength of selection for Party 1 relative to Party 2', type = float)
args = ap.parse_args()

# get arguments
if args.samples:
	samples = args.samples
else:
	samples = 500
# Define initial position and optima
position1 = np.zeros(1)
position = position1
position2 = position1
# Get divergence factor (default is 1)
if args.init_fit:
	r = 1-args.init_fit
else:
	r = 1-0.2
# Set average norm size for mutations
if args.ave_mut:
	average_mutation1 = args.ave_mut
	average_mutation2 = average_mutation1
	if args.rel_ave_mut:
		average_mutation1 = args.rel_ave_mut*average_mutation1
else:
	average_mutation1 = 0.1
	average_mutation2 = 0.1
# Get population sizes
# Population 1
if args.population_size1:
	N_1 = 10**(args.population_size1)
else:
	N_1 = 'infinite'
# Population 2
if args.population_size2:
	N_2 = 10**(args.population_size2)
else:
	N_2 = 'infinite'
# Get distributions
# Mutation distribution (default is uniform)
if args.mutations:
	distribution = args.mutations
else:
	distribution = 'normal'
# Number of mutations
if args.attempts:
	m = args.attempts
else:
	m = 5000
# Get mutation rate
if args.rate:
	rate = args.rate
else:
	rate = 1
if args.burn_in:
	burn_in = args.burn_in
else:
	burn_in = 0
if args.Q:
	Q = args.Q
else:
	Q = 2
if args.selection:
	d1 = args.selection
	d2 = d1
	if args.rel_selection:
		d2 = (1/args.rel_selection)*d2
else:
	d1 = 0.5
	d2 = 0.5

### OPEN OUTPUT FILES
output = open('FGM_data.csv', 'w')
output.write('Iteration,Simulation,z,s,Mutation Size,Fitness,Type,Population,Status\n')
# Ouputs mutations for FGMabiotic.py
mut_out = open('mut.csv', 'w')
mut_out.write('Mutation,Population\n')

### RUN SIMULATIONS
run_simulation(position, samples)
