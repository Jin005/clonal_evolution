'''

functions for non-spatial CSC model

author: Jacob G Scott 

21 Sep 2016

'''

import numpy as np

"""choose and define parent cell state"""
def choose_parent(state, stemflag):
	child_position = np.count_nonzero(state)
	parent_position = np.random.randint(child_position) #choose cell to divide
	parent = state[parent_position].astype(int)
	parentflag = stemflag[parent_position]
	return child_position, parent_position, parent, parentflag

def define_mutpair(mut_pairs, mut_number, num_muts, parent):
	parentname = parent.astype(str)
	child = (mut_number + num_muts).astype(int)
	childname = child[0].astype(str)
	mut_pairs.append((parentname, childname))
	mut_number = mut_number+num_muts
	return mut_pairs, mut_number, num_muts

"""discrete population dynamics, stem, non-spatial"""
def update_non_spatial_stem_CA(state, TACAge, TACAgeMax, stemflag, mut_pairs, mut_number, mut_rate, SymDiv):
	child_position, parent_position, parent, parentflag = choose_parent(state, stemflag)
	if parentflag == True: #if chosen cell is a stem cell
		if np.random.random() < SymDiv: #is this symmetric division?
			num_muts = np.random.poisson(mut_rate, 1).astype(int) #mutate?
			if num_muts == 0:
				state[child_position] = parent # new stem daughter (Tacked on end), is unmutated
				stemflag[child_position] = True
			else:
				mut_pairs, mut_number, num_muts = define_mutpair(mut_pairs, mut_number, num_muts, parent)
				state[child_position] = mut_number #add mutated cell at the end
				stemflag[child_position] = True
		else: #asymmetric division 	
			state[child_position] = parent # new TAC (on end)
			stemflag[child_position] = False
			num_muts = np.random.poisson(mut_rate, 1).astype(int) #get some mutations
			if num_muts > 0:
				mut_pairs, mut_number, num_muts = define_mutpair(mut_pairs, mut_number, num_muts, parent)
				state[parent_position] = mut_number #mutate parental cell
	else:#symmetric division of TAC
		state[child_position] = parent # new TAC (on end)
		stemflag[child_position] = False #new TAC
		TACAge[parent_position] += 1 #age parental TAC
		TACAge[child_position] = TACAge[parent_position] #child inherits parent age
		if TACAge[parent_position] == TACAgeMax:
			state = np.delete(state, parent_position)
			state = np.append(state, 0)
			stemflag = np.delete(stemflag, parent_position)
			stemflag = np.append(stemflag, False)
			TACAge = np.delete(TACAge, parent_position)
			TACAge = np.append(TACAge, 0)
	return state, TACAge, stemflag, mut_pairs, mut_number

def define_mutpair_list(state, mut_pairs, mut_number, num_muts, parent_mutation):
	mut_pairs.append((str(parent_mutation), str(mut_number + num_muts)))
	mut_number = mut_number+num_muts
	return mut_pairs, mut_number

"""discrete population dynamics, stem, non-spatial"""
def update_listbased_non_spatial_stem_CA(state, TACAgeMax, mut_pairs, mut_number, mut_rate, SymDiv):
	parent = np.random.randint(len(state))
	if state[parent][1] == -1: #if chosen cell is a stem cell
		if np.random.random() < SymDiv: #is this symmetric division?
			num_muts = np.random.poisson(mut_rate, 1)[0] #mutate?
			if num_muts == 0:
				state.append(state[parent]) # new stem daughter (Tacked on end), is unmutated
			else:
				mut_pairs, mut_number = define_mutpair_list(state, mut_pairs, mut_number, num_muts, state[parent][0])
				state.append((mut_number, -1)) #add mutated stem cell
		else: #asymmetric division 	
			state.append((state[parent][0], 0)) #add brand new TAC with parental mutation_number
			num_muts = np.random.poisson(mut_rate, 1)[0] #get some mutations
			if num_muts > 0:
				mut_pairs, mut_number = define_mutpair_list(state, mut_pairs, mut_number, num_muts, state[parent][0])
				state[parent] = (mut_number, -1) #mutate parental cell
	else:#symmetric division of TAC
		if state[parent][1] == TACAgeMax:
			state.pop(parent)
		else:
			state.append((state[parent][0], state[parent][1]+1)) #add new TAC with parental mutation_number and TACAge +1
			state[parent] = (state[parent][0], state[parent][1]+1) #modify parental TACAge
	return state, mut_pairs, mut_number

