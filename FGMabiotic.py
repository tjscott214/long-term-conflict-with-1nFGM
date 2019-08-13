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

# Function to generate random vectors

def generate_random_vector(n):
    # If simulation is 1d, randomly make radial positive or negative
    # Save vector as radial
    if n == 1:
    	# Define radial (magnitude of vector)
    	# Input determines which function to draw radial from
    	if distribution == 'uniform':
    		radial = np.random.uniform(0,uni)
    	if distribution == 'chi':
    		radial = np.random.chisquare(n)
    	if distribution == 'exponential':
    		radial = np.random.exponential(expo)
    	if distribution == 'normal':
    		radial = abs(np.random.normal(0, sd_1d))
    	vector = np.array(radial * (-1)**random.randint(1,2))
    # If simulation is > 1d, use spherical coordinates to construct vector
    if n > 1:
    	if distribution == 'normal':
    		mean = np.zeros(n)
    		sigma =  (np.identity(n)) * sd**2
    		vector = np.random.multivariate_normal(mean, sigma)
    		radial = dist.euclidean(mean, vector)
    	else:
    		# Define radial (magnitude of vector)
    		# Input determines which function to draw radial from
    		if distribution == 'uniform':
    			radial = np.random.uniform(0,10)
    		if distribution == 'chi':
    			# This distribution should mirror that of mutations drawn from a multivariate gaussian distribution
    			radial = np.random.chisquare(n)
    		if distribution == 'exponential':
    			radial = np.random.exponential(5)
    		# Create empty array for thetas
    		theta_array = np.array([])
    		# Create empty array for vector
    		vector = np.array([])
    		# Create product equal to 1
    		product = 1
    		# Initialize for loop over range n - 2
    		for i in range(n - 1):
    			# Calculate a random theta
        		theta = random.uniform(0, 2*np.pi)
        		# Append theta to array
        		theta_array = np.append(theta_array, theta)
    		# Calculate n - 1 theta randomly between 0 and pi and append to array
    		### Uncomment these to get unidirectional movement in two dimensions
    		#theta_n_minus_1 = np.pi * random.random()
    		#theta_array = np.append(theta_array, theta_n_minus_1)
 			# Loop over range n - 2
    		for i in range(n-2):
        		# Calculate cosine value for current iteration
        		new_cos = np.cos(theta_array[i])
        		# Calculate i component of vector
        		z = radial * product * new_cos
        		# Append to vector
        		vector = np.append(vector, z)
        		# Update product with sine of current theta
        		product = product * np.sin(theta_array[i])
    		# Calculate last two vector components with last theta in theta array
    		z_n_minus_1 = radial * product * np.cos(theta_array[-1])
    		vector = np.append(vector, z_n_minus_1)
    		z_n = radial * product * np.sin(theta_array[-1])
    		vector = np.append(vector, z_n)
    return radial, vector

def generate_optima(n1,n2):
	# Define the radial distance
	radial = (-2*np.log(r))**(1/Q)
	if args.top_one_axis:
		optimum1 = np.append(1, np.zeros(n1-1))
		optimum2 = np.append(-1, np.zeros(n2-1))
		jk = 'j'
	elif args.bottom_one_axis:
		optimum1 = np.append(np.zeros(n1-1), 1)
		optimum2 = np.append(np.zeros(n2-1), 1)
		jk = 'k'
	else:
		jk = '_'
		# Calculate z
		z1 = radial / (n1**(1/2)) 
		z2 = radial / (n2**(1/2)) 
		# Define an empty vector for the optimum
		optimum1 = np.array([])
		optimum2 = np.array([])
		# Loop through n and add z to vector
		for i in range(n1):
			optimum1 = np.append(optimum1, z1)
		for i in range(n2):
			optimum2 = np.append(optimum2, -z2)
		# Make optimum 2 the negative vector of optimum 1
	optimum2 = r2 * optimum2
	return optimum1, optimum2
# Function to norm distance
def norm_distance(distance):
	# divides distance by the standard deviation of a chi-square distribution with n degrees of freedom ()
	return distance
# Gaussian fitness function
def fitness_function(distance):
	return np.exp(-(((distance)**Q)/2))
# Define function that calculates probability of fixation event
def calculate_u(new_distance, old_distance, N = 'infinite'):
	# Calculate fitness of new mutation (e^(-d^2 / 2))
	fitness_new = fitness_function(new_distance)
	#print(fitness_new)
	# Calculate fitness of the old position
	fitness_old = fitness_function(old_distance)
	# Calculate s as relative fitness ((new fitness - old fitness) / old fitness)
	s_coefficient = (fitness_new/fitness_old) - 1
	# Use (1 - e^-2s)/(1 - e^-4Ns)
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
# Null model function
def shake(position, optimum, mut_list, order, generation_out, divergence_out):
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
	shake_optimum = optimum
	# Uncomment to examine per m mutations 
	for d in range(0,len(mut_list)):
	#while len(mutations) < 5:
		if order == 'Random':
			shake_optimum = shake_optimum + mut_list[random.randint(0,len(mut_list)-1)]*((-1)**(random.randint(1,2)))
		else:
			shake_optimum = shake_optimum + (mut_list[d])*((-1)**(random.randint(1,2)))
		#print(order, shake_optimum)
		distance_to_optimum = dist.euclidean(position, shake_optimum)
		# Create random vector
		mutation_size, vector = generate_random_vector(n)
		# Create new position
		future_position = position + vector
		# Calculate distance to optimum
		new_dist_to_optimum = dist.euclidean(future_position, shake_optimum)
		# Test whether distance is smaller than the probablility of fixation
		u, s = calculate_u(new_dist_to_optimum, distance_to_optimum, N_1)
		if random.random() <= u:
			mutation_fitness = vector
			# Update index
			divergence += 1
			# Update position
			position = future_position
			# Update distance to optimum
			distance_to_optimum = dist.euclidean(position, shake_optimum)
			if counter >= burn_in:
				fixed_out.write(str(s) + ',' + str(mutation_size) + ',Shake,' + order + ',Fixed\n')
		else:
			if counter >= burn_in:
				fixed_out.write(str(s) + ',' + str(mutation_size) + ',Shake,' + order + ',Unfixed\n')
		
		generational_fitness = distance_to_optimum
		generation_out.write(str(generational_fitness) + ',')
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
		# Change v is chase red queen is true
		if args.chase_red_queen:
			if random.random() <= switch_prob:
				v = random.randint(1,n-1)
		#print(v)
		# Create random vector from random vector function
		if random.randint(1,2) == 1:
			for j in range(rate):
				mutation_size1, random_vector = generate_random_vector(n)
				vector1 = random_vector
				# Add vector to position 
				future_position1 = position1 + random_vector
				# If n is 1, update population 2
				if n == 1:
					update_vector = random_vector
					future_position2 = future_position1
				# If n is greater than 1, update v number of axes
				elif args.chase_red_queen:
					#print(np.zeros(v-1), random_vector[v], np.zeros(n2-v))
					if v > 1:
						update_vector = np.append(np.zeros(v-1), np.array(random_vector[v]))
						update_vector = np.append(update_vector, np.zeros(n2 - v))
					else:
						update_vector = np.append(np.array(random_vector[v]), np.zeros(n2 - v))
					future_position2 = position2 + update_vector
				else:
					update_vector = np.append(random_vector[:v], np.zeros(n2 - v))
					future_position2 = position2 + update_vector
				#print(update_vector)
				# Calculate distances to optima
				new_dist_to_optimum1 = dist.euclidean(future_position1, optimum1)
				# Test whether distance is smaller than the probablility of fixation
				u1, s1 = calculate_u(new_dist_to_optimum1, distance_to_optimum1, N_1)
				if random.random() <= u1:
					# Add divergence
					divergence1 += 1
					position1 = future_position1
					position2 = future_position2
					# Update mutation fitness
					mutation_fitness1_1 = random_vector
					mutation_fitness2_1 = update_vector
					#mutation_fitness1_1 = fitness_function(dist.euclidean(position1, optimum1))/w0_1
					#mutation_fitness2_1 = fitness_function(dist.euclidean(position2, optimum2))/w0_2
					# Update distance to optimum
					distance_to_optimum1 = new_dist_to_optimum1
					distance_to_optimum2 = dist.euclidean(position2, optimum2)
					# Update mutation fitness
					#mutation_fitness1_1 = fitness_function(distance_to_optimum1)
					#mutation_fitness2_1 = fitness_function(distance_to_optimum2)
					if counter >= burn_in:
						fixed_out.write(str(s1) + ',' + str(mutation_size1) + ',Conflict,Population 1,Fixed,' + str(v) + '\n')
				else:
					if counter >= burn_in:
						fixed_out.write(str(s1) + ',' + str(mutation_size1) + ',Conflict,Population 1,Unfixed,' + str(v) + '\n')
			### Population 2 mutation
			mutation_size2, random_vector = generate_random_vector(n2)
			# Add vector to position 
			future_position2 = position2 + random_vector
			vector2 = random_vector
			# If n is 1, update population 1
			if n == 1:
				update_vector = random_vector
				future_position1 = future_position2
			# If n is greater than 1, update v number of axes
			elif args.chase_red_queen:
				#print(np.zeros(v-1), random_vector[v], np.zeros(n-v))
				if v > 1:
					update_vector = np.append(np.zeros(v-1), np.array(random_vector[v]))
					update_vector = np.append(update_vector, np.zeros(n - v))
				else:
					update_vector = np.append(np.array(random_vector[v]), np.zeros(n - v)) 
				future_position1 = position1 + update_vector
			else:
				update_vector = np.append(random_vector[:v], np.zeros(n - v))
				future_position1 = position1 + update_vector
			# Calculate distances to optima
			new_dist_to_optimum2 = dist.euclidean(future_position2, optimum2)
			# Test whether distance is smaller than the probablility of fixation
			u2, s2 = calculate_u(new_dist_to_optimum2, distance_to_optimum2, N_2)
			if random.random() <= u2:
				divergence2 += 1
				position1 = future_position1
				position2 = future_position2
				# Update mutation fitness
				mutation_fitness1_2 = update_vector
				mutation_fitness2_2 = random_vector
				#mutation_fitness1_2 = fitness_function(dist.euclidean(position1, optimum1))/w0_1
				#mutation_fitness2_2 = fitness_function(dist.euclidean(position2, optimum2))/w0_2
				# Update distance to optimum
				distance_to_optimum1 = dist.euclidean(position1, optimum1)
				distance_to_optimum2 = new_dist_to_optimum2
				# Update mutation fitness
				#mutation_fitness1_2 = fitness_function(distance_to_optimum1)
				#mutation_fitness2_2 = fitness_function(distance_to_optimum2)
				if counter >= burn_in:
					fixed_out.write(str(s2) + ',' + str(mutation_size2) + ',Conflict,Population 2,Fixed,' + str(v) + '\n')
			else:
				if counter >= burn_in:
					fixed_out.write(str(s2) + ',' + str(mutation_size2) + ',Conflict,Population 2,Unfixed,' + str(v) + '\n')
			counter += 1
			if n == 1:
				generation_out1.write(str(position1[0]) + ',')
				generation_out2.write(str(position2[0]) + ',')
			else:	
				generational_fitness1 = distance_to_optimum1
				generational_fitness2 = distance_to_optimum2
				generation_out1.write(str(generational_fitness1) + ',')
				generation_out2.write(str(generational_fitness2) + ',')
			divergence_out1.write(str(divergence1) + ',')
			divergence_out2.write(str(divergence2) + ',')
		else:
			### Population 2 mutates
			mutation_size2, random_vector = generate_random_vector(n2)	
			vector2 = random_vector
			# Add vector to position 
			future_position2 = position2 + random_vector
			# If n is 1, update population 1
			if n == 1:
				update_vector = random_vector
				future_position1 = future_position2
			# If n is greater than 1, update v number of axes
			elif args.chase_red_queen:
				#print(np.zeros(v-1), random_vector[v], np.zeros(n-v))
				if v > 1:
					update_vector = np.append(np.zeros(v-1), np.array(random_vector[v]))
					update_vector = np.append(update_vector, np.zeros(n - v))
				else:
					update_vector = np.append(np.array(random_vector[v]), np.zeros(n - v))
				future_position1 = position1 + update_vector
			else:
				update_vector = np.append(random_vector[:v], np.zeros(n - v))
				future_position1 = position1 + update_vector
			# Calculate distances to optima
			new_dist_to_optimum2 = dist.euclidean(future_position2, optimum2)
			# Test whether distance is smaller than the probablility of fixation
			u2, s2 = calculate_u(new_dist_to_optimum2, distance_to_optimum2, N_2)
			if random.random() <= u2:
				divergence2 += 1
				position1 = future_position1
				position2 = future_position2
				# Update mutation fitness
				mutation_fitness1_1 = update_vector
				mutation_fitness2_1 = random_vector
				#mutation_fitness1_2 = fitness_function(dist.euclidean(position1, optimum1))/fitness_function(distance_to_optimum1)
				#mutation_fitness2_1 = fitness_function(dist.euclidean(position2, optimum2))/fitness_function(distance_to_optimum2)
				# Update distance to optimum
				distance_to_optimum1 = dist.euclidean(position1, optimum1)
				distance_to_optimum2 = new_dist_to_optimum2
				# Update mutation fitness
				#mutation_fitness1_1 = fitness_function(distance_to_optimum1)
				#mutation_fitness2_1 = fitness_function(distance_to_optimum2)
				if counter >= burn_in:
					fixed_out.write(str(s2) + ',' + str(mutation_size2) + ',Conflict,Population 2,Fixed,' + str(v) + '\n')
			else:
				if counter >= burn_in:
					fixed_out.write(str(s2) + ',' + str(mutation_size2) + ',Conflict,Population 2,Unfixed,' + str(v) + '\n')
			### Population 1 mutates
			for j in range(rate):
				mutation_size1, random_vector = generate_random_vector(n)
				vector1 = random_vector
				# Add vector to position 
				future_position1 = position1 + random_vector
				# If n is 1, update population 2
				if n == 1:
					update_vector = random_vector
					future_position2 = future_position1
				# If n is greater than 1, update v number of axes
				elif args.chase_red_queen:
					#print(np.zeros(v-1), random_vector[v], np.zeros(n2-v))
					if v > 1:
						update_vector = np.append(np.zeros(v-1), np.array(random_vector[v]))
						update_vector = np.append(update_vector, np.zeros(n2 - v))
					else:
						update_vector = np.append(np.array(random_vector[v]), np.zeros(n2 - v))
					future_position2 = position2 + update_vector
				else:
					update_vector = np.append(random_vector[:v], np.zeros(n2 - v))
					future_position2 = position2 + update_vector
				# Calculate distances to optima
				new_dist_to_optimum1 = dist.euclidean(future_position1, optimum1)
				# Test whether distance is smaller than the probablility of fixation
				u1, s1 = calculate_u(new_dist_to_optimum1, distance_to_optimum1, N_1)
				if random.random() <= u1:
					divergence1 += 1
					position1 = future_position1
					position2 = future_position2
					# Update mutation fitness
					mutation_fitness1_2 = random_vector
					mutation_fitness2_2 = update_vector
					#mutation_fitness1_1 = fitness_function(dist.euclidean(position1, optimum1))/fitness_function(distance_to_optimum1)
					#mutation_fitness2_2 = fitness_function(dist.euclidean(position2, optimum2))/fitness_function(distance_to_optimum2)
					# Update distance to optimum
					distance_to_optimum1 = new_dist_to_optimum1
					distance_to_optimum2 = dist.euclidean(position2, optimum2)
					# Update mutation fitness
					#mutation_fitness1_2 = fitness_function(distance_to_optimum1)
					#mutation_fitness2_2 = fitness_function(distance_to_optimum2)
					if counter >= burn_in:
						fixed_out.write(str(s1) + ',' + str(mutation_size1) + ',Conflict,Population 1,Fixed,' + str(v) + '\n')
				else:
					if counter >= burn_in:
						fixed_out.write(str(s1) + ',' + str(mutation_size1) + ',Conflict,Population 1,Unfixed,' + str(v) + '\n')
			counter += 1
			if n == 1:
				generation_out1.write(str(position1[0]) + ',')
				generation_out2.write(str(position2[0]) + ',')
			else:	
				generational_fitness1 = distance_to_optimum1
				generational_fitness2 = distance_to_optimum2
				generation_out1.write(str(generational_fitness1) + ',')
				generation_out2.write(str(generational_fitness2) + ',')
			divergence_out1.write(str(divergence1) + ',')
			divergence_out2.write(str(divergence2) + ',')
	# Add newlines
	generation_out1.write('\n')
	generation_out2.write('\n')
	divergence_out1.write('\n')
	divergence_out2.write('\n')
def resample(position, num_samples,v):
	# Create a numpy array for fixations
	fixations_array = np.array([])
	fixations1_array = np.array([])
	fixations2_array = np.array([])
	df = pd.read_csv(shake_file)
	optimum1, optimum2 = generate_optima(n,n2)
	#master_mut_list = df['Mutation'].tolist()
	master_mut_list = df.groupby('Population')['Mutation'].apply(list)[1]
	#print(master_mut_list[0:200])
	max_0 = max(sum(1 for _ in g) for k, g in groupby(master_mut_list) if k==0)
	print("Proportion of zeros: ", master_mut_list.count(0)/len(master_mut_list), 'Max cons: ', max_0)
	#print(optimum1, optimum2)
	# Resample using for loop
	index = 0 
	for sample in range(num_samples):
		mut_list = master_mut_list[index:index + m]
		#print(index, index+m, m)
		shake(position, optimum1, mut_list, 'Shake', generation_out_ordered, divergence_out_ordered)
		shake(position, optimum1, mut_list, 'Random', generation_out_random, divergence_out_random)
		index += m
	# Close output files
	generation_out_ordered.close()
	generation_out_random.close()
	#generation_out1.close()
	#generation_out2.close()
	divergence_out_ordered.close()
	divergence_out_random.close()
	#divergence_out1.close()
	#divergence_out2.close()
	fixed_out.close()





# Parse arguments using argparse
ap = argparse.ArgumentParser()
ap.add_argument('-x', '--samples', help = 'number of resamples', type = int)
ap.add_argument('-p', '--population_size1', help = 'population size for one population', type = int)
ap.add_argument('-pp', '--population_size2', help = 'population size for second population', type = int)
ap.add_argument('-m', '--mutations', help = 'mutation distribution for mutation vectors')
ap.add_argument('-q', '--Q', help = 'changes Q parameter in fitness function', type = float)
ap.add_argument('-z', '--attempts', help = 'number of generations per walk', type = int)
ap.add_argument('-c', '--init_fit', help = 'changes the distance optimal values by a factor of the input value', type = float)
ap.add_argument('-cc', '--factor2', help = 'changes the distance of the negative optimal value by a factor of the input value', type = float)
ap.add_argument('-r', '--rate', help = 'mutation rate for population 1', type = int)
ap.add_argument('-b', '--burn_in', help = 'define burn in period for equilibrium', type = int)
ap.add_argument('-a', '--ave_mut', help = 'average mutation norm', type = float)
ap.add_argument('-mut', '--shake', help = 'mutation file for shaking optimum', type = str)
args = ap.parse_args()

# Define n, number of joint traits, and samples
n = args.traits
n_base = str(args.traits) + 'n'
if args.alt_traits:
	n2 = args.alt_traits
	n2_base = str(args.alt_traits) + 'n2'
else:
	n2 = n
	n2_base = ''
if args.joint_traits:
	v = args.joint_traits
else:
	v = 1
v_base =  '_v' + str(v) + '_'
if args.samples:
	samples = args.samples
else:
	samples = 500
# Normalization factor that affects fitness function and position of optimum
norm_factor = dist.euclidean(np.ones(n), np.zeros(n))
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
else:
	r = 1
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
ax = 10
# Calculate normalization factor
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
if args.shake:
	shake_file = args.shake[:-7] + 'mut.csv'
base_file = q_string + str(N_1) + '_' + str(N_2) + '_' + r_strng + 'df_' + str(r2) + 'ms_' + average_mutation_string + '_' + '_x' + str(rate) + jk + str(burn_in) + 'burn_' + crq_string
# Open files for generational data
generation_out_ordered = open(base_file + 'shake_gen.csv', 'w')
generation_out_random = open(base_file + 'random_gen.csv', 'w')	
#generation_out1 = open(base_file + 'gen1.csv', 'w')	
#generation_out2 = open(base_file + 'gen2.csv', 'w')
divergence_out_ordered = open(base_file + 'shake_div.csv', 'w')	
divergence_out_random = open(base_file + 'random_div.csv', 'w')	
#divergence_out1 = open(base_file + 'div1.csv', 'w')	
#divergence_out2 = open(base_file + 'div2.csv', 'w')	
fixed_out = open(base_file + 'shake_data.csv', 'w')
fixed_out.write('s,Mutation Size,Type,Population,Status,v\n')
#unfixed_out = open(n_base + v_base + str(N_1) + '_' + str(N_2) + '_'+ str(r) + 'df_' + str(r2) + 'ms_' + str(average_mutation) + '_' + '_x' + str(rate) + jk + str(burn_in) + 'burn_' + 'unfixed.csv', 'w')
#unfixed_out.write('s,Mutation Size,Type,Population\n')
### Run simulations
resample(position, samples, v)
