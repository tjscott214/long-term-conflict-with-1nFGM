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
# Function to generate random vectors
def generate_random_vector():
    if distribution == 'uniform':
    	radial = np.random.uniform(0,uni)
    if distribution == 'exponential':
    	radial = np.random.exponential(expo)
    if distribution == 'normal':
    	radial = abs(np.random.normal(0, sd_1d))
    vector = np.array(radial * (-1)**random.randint(1,2))
    # Return radial and vector
    return radial, vector

def generate_optima(d1, d2):
	# Define the intersection distance
	optimum1 = np.array([(-d1*np.log(r))**(1/Q)])
	optimum2 = np.array([-(-d2*np.log(r))**(1/Q)])
	# Return optimum values
	return optimum1, optimum2
# Gaussian fitness function
def fitness_function(distance, denominator):
	return np.exp(-(((distance)**Q)/denominator))
# Define function that calculates probability of fixation event
def calculate_u(new_distance, old_distance, N = 'infinite', denominator = 2):
	# Calculate fitness of new mutation (e^(-d^2 / 2))
	fitness_new = fitness_function(new_distance, denominator)
	# Calculate fitness of the old position
	fitness_old = fitness_function(old_distance, denominator)
	# Calculate s as relative fitness ((new fitness - old fitness) / old fitness)
	s_coefficient = (fitness_new/fitness_old) - 1
	# Use (1 - e^-2s) or (1 - e^-2s)/(1 - e^-4Ns)
	# Calculate probablility of fixation
	if N == 'infinite':
		probability_of_fixation = (1 - np.exp(-2*s_coefficient))
	elif N > 0:
		probability_of_fixation = ((1 - np.exp(-2*s_coefficient))/(1 - np.exp(-4*s_coefficient*N)))
	# Return the probability that a new mutation will be fixed
	return probability_of_fixation, s_coefficient
# Calculate individual distances
def calculate_axes(position, optimum):
	# Distances array
	distances = np.array([])
	# Loop through position indices
	for i in range(len(position)):
		distance = dist.euclidean(position[i], optimum[i])
		# Append distance to distances array
		distances = np.append(distances, distance)
	#return array
	return distances
# Non conflict model function
def non_conflict(position, optimum):
	# Define mutation counter
	counter = 0
	#Define divergence counter
	divergence = 0
	# Define distance to optimum
	distance_to_optimum = dist.euclidean(position, optimum)
	# Create empty list for dataframes
	distances_list = []
	fitness_list = []
	position_list = []
	# Simulate for m iterations
	while counter < m:
	#while len(mutations) < 5:
		# Create random vector
		mutation_size, vector = generate_random_vector()
		#distribution = np.append(distribution, mutation_size)
		# Create new position
		future_position = position + vector
		# Calculate distance to optimum
		new_dist_to_optimum = dist.euclidean(future_position, optimum)
		# Test whether distance is smaller than the probablility of fixation
		u, s = calculate_u(new_dist_to_optimum, distance_to_optimum, N_1, d1)
		if random.random() <= u:
			mutation_fitness = vector
			# Update index
			divergence += 1
			# Update position
			position = future_position
			# Update distance to optimum
			distance_to_optimum = new_dist_to_optimum
			if counter >= burn_in:
				fixed_out.write(str(s) + ',' + str(mutation_size) + ',No Conflict,No Conflict,Fixed\n')
		else:
			if counter >= burn_in:
				fixed_out.write(str(s) + ',' + str(mutation_size) + ',No Conflict,No Conflict,Unfixed\n')
		generation_out.write(str(position[0]) + ',')
		divergence_out.write(str(divergence)+ ',')
		counter += 1
	# Add new line to generation file
	generation_out.write('\n')
	divergence_out.write('\n')
# Two population function
def two_population_model(position1, position2, optimum1, optimum2,v):
	# Define starting positions
	position1 = position1
	position2 = position2 
	# Define mutation counter
	counter = 0
	# Define divergence
	divergence1 = 0
	divergence2 = 0
	# Define distance to optima 
	# Distance is calculated using a scipy distance function
	distance_to_optimum1 = dist.euclidean(position1, optimum1)
	distance_to_optimum2 = dist.euclidean(position2, optimum2)
	# Create empty list for dataframes
	# Create empty list for dataframes
	distances1_list = []
	distances2_list = []
	fitness1_list = []
	fitness2_list = []
	position1_list = []
	position2_list = []
	#print(v)
	# Loop until mutations exceed m
	while counter < m:
		# Test which party mutates first
		if random.randint(1,2) == 1:
			# Party 1 mutates first
			for j in range(rate):
				#print(j)
				mutation_size1, random_vector = generate_random_vector()
				vector1 = random_vector
				# Add vector to position 
				future_position1 = position1 + random_vector
				# If n is 1, update population 2
				update_vector = random_vector
				future_position2 = future_position1
				# Calculate distances to optima
				new_dist_to_optimum1 = dist.euclidean(future_position1, optimum1)
				# Test whether distance is smaller than the probablility of fixation
				u1, s1 = calculate_u(new_dist_to_optimum1, distance_to_optimum1, N_1, d1)
				if random.random() <= u1:
					# Add divergence
					divergence1 += 1
					position1 = future_position1
					position2 = future_position2
					# Update mutation fitness
					mutation_fitness1_1 = random_vector
					mutation_fitness2_1 = update_vector
					# Update distance to optimum
					distance_to_optimum1 = new_dist_to_optimum1
					distance_to_optimum2 = dist.euclidean(position2, optimum2)
					# Output fixation data
					mut_out.write(str(vector1) + ',Population 1\n')
					if counter >= burn_in:
						fixed_out.write(str(s1) + ',' + str(mutation_size1) + ',Conflict,Population 1,Fixed,' + str(v) + '\n')
				else:
					# Output failed fixation data
					if counter >= burn_in:
						fixed_out.write(str(s1) + ',' + str(mutation_size1) + ',Conflict,Population 1,Unfixed,' + str(v) + '\n')
					mut_out.write('0,Population 1\n')
			### Party 2 mutation
			mutation_size2, random_vector = generate_random_vector()
			# Add vector to position 
			future_position2 = position2 + random_vector
			vector2 = random_vector
			# If n is 1, update population 1
			update_vector = random_vector
			future_position1 = future_position2
			# Calculate distances to optima
			new_dist_to_optimum2 = dist.euclidean(future_position2, optimum2)
			# Test whether distance is smaller than the probablility of fixation
			u2, s2 = calculate_u(new_dist_to_optimum2, distance_to_optimum2, N_2, d2)
			if random.random() <= u2:
				divergence2 += 1
				position1 = future_position1
				position2 = future_position2
				# Update mutation fitness
				mutation_fitness1_2 = update_vector
				mutation_fitness2_2 = random_vector
				# Update distance to optimum
				distance_to_optimum1 = dist.euclidean(position1, optimum1)
				distance_to_optimum2 = new_dist_to_optimum2
				# Output fixation data
				mut_out.write(str(vector2) + ',Population 2\n')
				if counter >= burn_in:
					fixed_out.write(str(s2) + ',' + str(mutation_size2) + ',Conflict,Population 2,Fixed,' + str(v) + '\n')
			else:
				# Output failed fixation data
				if counter >= burn_in:
					fixed_out.write(str(s2) + ',' + str(mutation_size2) + ',Conflict,Population 2,Unfixed,' + str(v) + '\n')
				mut_out.write('0,Population 2\n')
			counter += 1
			generation_out1.write(str(position1[0]) + ',')
			generation_out2.write(str(position2[0]) + ',')
			divergence_out1.write(str(divergence1) + ',')
			divergence_out2.write(str(divergence2) + ',')
		# Scenario where party 2 starts
		else:
			### Party 2 mutates
			mutation_size2, random_vector = generate_random_vector()	
			vector2 = random_vector
			# Add vector to position 
			future_position2 = position2 + random_vector
			# If n is 1, update population 1
			update_vector = random_vector
			future_position1 = future_position2
			# Calculate distances to optima
			new_dist_to_optimum2 = dist.euclidean(future_position2, optimum2)
			# Test whether distance is smaller than the probablility of fixation
			u2, s2 = calculate_u(new_dist_to_optimum2, distance_to_optimum2, N_2, d2)
			if random.random() <= u2:
				divergence2 += 1
				position1 = future_position1
				position2 = future_position2
				# Update mutation fitness
				mutation_fitness1_1 = update_vector
				mutation_fitness2_1 = random_vector
				# Update distance to optimum
				distance_to_optimum1 = dist.euclidean(position1, optimum1)
				distance_to_optimum2 = new_dist_to_optimum2
				# Output fixation data
				mut_out.write(str(vector2) + ',Population 2\n')
				if counter >= burn_in:
					fixed_out.write(str(s2) + ',' + str(mutation_size2) + ',Conflict,Population 2,Fixed,' + str(v) + '\n')
			else:
				# Output failed fixation data
				if counter >= burn_in:
					fixed_out.write(str(s2) + ',' + str(mutation_size2) + ',Conflict,Population 2,Unfixed,' + str(v) + '\n')
				mut_out.write('0,Population 2\n')
			### Party 1 mutates
			for j in range(rate):
				mutation_size1, random_vector = generate_random_vector()
				vector1 = random_vector
				# Add vector to position 
				future_position1 = position1 + random_vector
				# If n is 1, update population 2
				update_vector = random_vector
				future_position2 = future_position1
				# Calculate distances to optima
				new_dist_to_optimum1 = dist.euclidean(future_position1, optimum1)
				# Test whether distance is smaller than the probablility of fixation
				u1, s1 = calculate_u(new_dist_to_optimum1, distance_to_optimum1, N_1, d1)
				if random.random() <= u1:
					divergence1 += 1
					position1 = future_position1
					position2 = future_position2
					# Update mutation fitness
					mutation_fitness1_2 = random_vector
					mutation_fitness2_2 = update_vector
					# Update distance to optimum
					distance_to_optimum1 = new_dist_to_optimum1
					distance_to_optimum2 = dist.euclidean(position2, optimum2)
					# Output fixation data
					mut_out.write(str(vector1) + ',Population 1\n')
					if counter >= burn_in:
						fixed_out.write(str(s1) + ',' + str(mutation_size1) + ',Conflict,Population 1,Fixed,' + str(v) + '\n')
				else:
					# Output failed fixation data
					if counter >= burn_in:
						fixed_out.write(str(s1) + ',' + str(mutation_size1) + ',Conflict,Population 1,Unfixed,' + str(v) + '\n')
					mut_out.write('0,Population 1\n')
			counter += 1
			generation_out1.write(str(position1[0]) + ',')
			generation_out2.write(str(position2[0]) + ',')
			divergence_out1.write(str(divergence1) + ',')
			divergence_out2.write(str(divergence2) + ',')
	# Add newlines
	generation_out1.write('\n')
	generation_out2.write('\n')
	divergence_out1.write('\n')
	divergence_out2.write('\n')
# Resamples single simulations
def run_simulation(position, num_samples,v):
	# Define optimum values
	optimum1, optimum2 = generate_optima(n,d1,d2)
	print(optimum1, optimum2)
	print(fitness_function(dist.euclidean(optimum1,0), d1), fitness_function(dist.euclidean(optimum2,0), d2))
	# Resample using for loop
	for sample in range(num_samples):
		non_conflict(position, optimum1)
		two_population_model(position1, position2, optimum1, optimum2,v)
	# Close output files
	generation_out.close()
	generation_out1.close()
	generation_out2.close()
	divergence_out.close()
	divergence_out1.close()
	divergence_out2.close()
	fixed_out.close()
	mut_out.close()
	#unfixed_out.close()

### SET ARGUMENTS
ap = argparse.ArgumentParser()
ap.add_argument('-x', '--samples', help = 'number of resamples', type = int)
ap.add_argument('-p', '--population_size1', help = 'population size for first party', type = int)
ap.add_argument('-pp', '--population_size2', help = 'population size for second party', type = int)
ap.add_argument('-m', '--mutations', help = 'mutation distribution for mutation vectors (normal-default, uniform, or exponential)')
ap.add_argument('-q', '--Q', help = 'changes Q (epistasis) parameter in fitness function', type = float)
ap.add_argument('-z', '--attempts', help = 'number of iterations per single simulation', type = int)
ap.add_argument('-c', '--init_fit', help = 'defines initial fitness', type = float)
ap.add_argument('-cc', '--factor2', help = 'multiplies only the second party\'s optimum', type = float)
ap.add_argument('-r', '--rate', help = 'allows first party to mutate r times for every one mutation of second party', type = int)
ap.add_argument('-b', '--burn_in', help = 'defines burn in period for equilibrium', type = int)
ap.add_argument('-a', '--ave_mut', help = 'scales average mutation sizes', type = float)
ap.add_argument('-d', '--selection', help = 'Adjust strength of selection for Party 1 relative to Party 2', type = float)
args = ap.parse_args()


if args.samples:
	samples = args.samples
else:
	samples = 500
# Define initial position and optima
position1 = np.zeros(n)
position = position1
position2 = np.zeros(n2)
# Get divergence factor (default is 1)
if args.init_fit:
	r = args.init_fit
	if len(str(r)) == 3:
		r_strng = str(r) + '0'
	else:
		r_strng = str(r)
if args.factor2:
	r2 = args.factor2
else:
	r2 = 1
# Set average norm size for mutations
if args.ave_mut:
	average_mutation = args.ave_mut
	if len(str(average_mutation)) == 3:
		average_mutation_string = str(average_mutation) + '0'
	else:
		average_mutation_string = str(average_mutation)
else:
	average_mutation = 0.5
	average_mutation_string = '0.50'
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
	print(rate)
else:
	rate = 1
ax = 10
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
	d1 = 2
	d2 = args.selection*2
else:
	d1 = 2
	d2 = 2
# Calculate scaling factors for normal, uniform, and exponential
sd_1d = average_mutation*((np.pi)**(1/2))/(2**(1/2))
uni = 2*average_mutation
expo = average_mutation

### OPEN OUTPUT FILES
# Base file name
base_file = q_string + str(N_1) + '_' + str(N_2) + '_' + str(r_strng) + 'df_' + str(r2) + 'ms_' + average_mutation_string + '_' + '_x' + str(rate) + '_' + str(burn_in) + 'burn_'
# Open files for output data
# Outputs columns of position values for each iteration
generation_out = open(base_file + 'gen.csv', 'w')	
generation_out1 = open(base_file + 'gen1.csv', 'w')	
generation_out2 = open(base_file + 'gen2.csv', 'w')
# Outputs columns of divergence (number of fixations) values for each iteration
divergence_out = open(base_file + 'div.csv', 'w')	
divergence_out1 = open(base_file + 'div1.csv', 'w')	
divergence_out2 = open(base_file + 'div2.csv', 'w')	
# Outputs selection coefficients, mutation sizes
fixed_out = open(base_file + 'data.csv', 'w')
fixed_out.write('s,Mutation Size,Type,Population,Status,v\n')
# Ouputs mutations by axis
mut_out = open(base_file + 'mut.csv', 'w')
if n == 1:
	mut_out.write('Mutation,Population\n')
else:
	for j in range(n - 1):
		mut_out.write('Trait ' + str(j) + ',')
	mut_out.write('Trait ' + str(n) + ',Population\n')

### RUN SIMULATIONS
run_simulation(position, samples, v)
