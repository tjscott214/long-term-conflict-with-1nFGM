#!/usr/bin/env python

### This program simulates Fisher's geometric model with abiotic change equal to fixations during conflict simulations (from FGMconflict.py) ###
### python3 FGMabiotic.py -help for input options ###
### Written by Trey J Scott 2018 ###
### python --version ###
### Python 3.5.2 :: Anaconda 4.2.0 (x86_64) ###

# Import programs
import random
import numpy as np
from scipy.spatial import distance as dist
from scipy.stats import norm
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import scipy.special as spc
from itertools import groupby

### FUNCTIONS ###

# Function to generate random mutations with a specified average size
def generate_random_vector():
    if distribution == 'uniform':
    	radial = np.random.uniform(0,uni)
    if distribution == 'chi':
    	radial = np.random.chisquare(n)
    if distribution == 'exponential':
    	radial = np.random.exponential(expo)
    if distribution == 'normal':
    	radial = abs(np.random.normal(0, sd_1d))
    vector = np.array(radial * (-1)**random.randint(1,2))
    return radial, vector

# Gaussian fitness function
def fitness_function(distance,d):
	return np.exp(-(d*(distance**Q)))

# Calculates probability of fixation for new mutations
def calculate_u(new_distance, old_distance, N = 'infinite', denominator = 0.5):
	fitness_new = fitness_function(new_distance, denominator)
	fitness_old = fitness_function(old_distance, denominator)
	s_coefficient = (fitness_new/fitness_old) - 1
	if N == 'infinite':
		probability_of_fixation = (1 - np.exp(-2*s_coefficient))
	elif N > 0:
		probability_of_fixation = ((1 - np.exp(-2*s_coefficient))/(1 - np.exp(-4*s_coefficient*N)))
	return probability_of_fixation, s_coefficient

# Functon that simulates adaptation to a moving optimum with Fisher's geometric model
def abiotic_change(position, optimum, mut_list, samp):
	counter = 0
	distance_to_optimum = dist.euclidean(position, optimum)
	moving_optimum = optimum
	for d in range(0,len(mut_list)):
		moving_optimum = moving_optimum + (mut_list[d])*((-1)**(random.randint(1,2)))
		distance_to_optimum = dist.euclidean(position, moving_optimum)
		mutation_size, vector = generate_random_vector()
		future_position = position + vector
		new_dist_to_optimum = dist.euclidean(future_position, moving_optimum)
		u, s = calculate_u(new_dist_to_optimum, distance_to_optimum, N_1,d1)
		if random.random() <= u:
			mutation_fitness = vector
			position = future_position
			distance_to_optimum = dist.euclidean(position, moving_optimum)
			if counter >= burn_in:
				output.write(str(d) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s) + ',' + str(mutation_size) + ',' + str(fitness_function(distance_to_optimum,d1)) + ',Abiotic Change,Fixed\n')
		else:
			if counter >= burn_in:
				output.write(str(d) + ',' + str(samp) + ',' + str(position[0]) + ',' + str(s) + ',' + str(mutation_size)+ ',' + str(fitness_function(distance_to_optimum,d1)) + ',Abiotic Change,Unfixed\n')
		counter += 1

# Runs simulations multiple times
def run_simulations(position, num_samples):
	df = pd.read_csv(shake_file)
	optimum = np.array([(-(1/d1)*np.log(r))**(1/Q)])
	master_mut_list = df.groupby('Population')['Mutation'].apply(list)[1]
	index = 0 
	for sample in range(num_samples):
		mut_list = master_mut_list[index:index + m]
		abiotic_change(position, optimum, mut_list, sample)
		index += m
	output.close()





### SET ARGUMENTS
ap = argparse.ArgumentParser()
ap.add_argument('-x', '--samples', help = 'number of resamples', type = int)
ap.add_argument('-p', '--population_size1', help = 'population size for one population', type = int)
ap.add_argument('-pp', '--population_size2', help = 'population size for second population', type = int)
ap.add_argument('-m', '--mutations', help = 'mutation distribution for mutation vectors')
ap.add_argument('-q', '--Q', help = 'changes Q parameter in fitness function', type = float)
ap.add_argument('-z', '--attempts', help = 'number of generations per walk', type = int)
ap.add_argument('-c', '--init_fit', help = 'changes the distance optimal values by a factor of the input value', type = float)
ap.add_argument('-r', '--rate', help = 'mutation rate for population 1', type = int)
ap.add_argument('-b', '--burn_in', help = 'define burn in period for equilibrium', type = int)
ap.add_argument('-a', '--ave_mut', help = 'average mutation norm', type = float)
ap.add_argument('-d', '--selection', help = 'Adjust strength of selection', type = float)
ap.add_argument('-mut', '--changes', help = 'mutation file for moving optimum', type = str)
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

if args.init_fit:
	r = 1-args.init_fit
else:
	r = 1-0.2
# Set average norm size for mutations
if args.ave_mut:
	average_mutation = args.ave_mut
else:
	average_mutation = 0.1
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
	m = 50000
# Get mutation rate
if args.rate:
	rate = args.rate
else:
	rate = 1
# Calculate normalization factor (used in mutation function)
sd_1d = average_mutation*((np.pi)**(1/2))/(2**(1/2))
uni = 2*average_mutation
expo = average_mutation
if args.burn_in:
	burn_in = args.burn_in
else:
	burn_in = 0
if args.Q:
	Q = args.Q
	q_string = 'Q_' + str(Q) + '_'
else:
	Q = 2
	q_string = ''
if args.selection:
	d1 = args.selection
else:
	d1 = 0.5
if args.changes:
	shake_file = args.changes[:-7] + 'mut.csv'

# Open output file
output = open('abiotic_data.csv', 'w')
output.write('Iteration,Simulation,z,s,Mutation Size,Fitness,Population,Status\n')

### Run simulations
run_simulations(position, samples)
