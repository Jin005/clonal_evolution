'''
non-spatial CA to determine effects of spatial organization on phylogenies

july 17th, 2016

jacob g. scott
'''
# import dendropy
# from dendropy import Tree as DTree
# from dendropy.calculate import treemeasure
import numpy as np
import matplotlib.pyplot as plt
from random import random
# from ete3 import Tree
import clonal_evolution_functions as cef
import pandas as pd
import seaborn as sns
from pylab import rcParams


"""discrete population dynamics, stem, non-spatial"""
def update_non_spatial_stem_CA(state, TACAge, stemflag, mut_pairs, mut_number, SymDiv):
	child_position = np.count_nonzero(state)
	parent_position = np.random.randint(child_position) #choose cell to divide
	parent = state[parent_position].astype(int)
	parentflag = stemflag[parent_position]
	if parentflag == True: #if chosen cell is a stem cell
		if random() < SymDiv: #is this symmetric division?
			num_muts = np.random.poisson(mut_rate, 1).astype(int) #get some mutations
			if num_muts == 0:
				state[child_position] = parent # new stem daughter (Tacked on end), is unmutated
				stemflag[child_position] = True
			else:
				parentname = parent.astype(str)
				child = (mut_number + num_muts).astype(int)
				childname = child[0].astype(str)
				mut_pairs.append((parentname, childname))
				mut_number = mut_number+num_muts
				state[child_position] = mut_number #add mutated cell at the end
				stemflag[child_position] = True
		else: #asymmetric division 	
			state[child_position] = parent # new TAC (on end)
			stemflag[child_position] = False
			num_muts = np.random.poisson(mut_rate, 1).astype(int) #get some mutations
			if num_muts > 0:
				parentname = parent.astype(str)
				child = (mut_number + num_muts).astype(int)
				childname = child[0].astype(str)
				mut_pairs.append((parentname, childname))
				mut_number = mut_number + num_muts
				state[parent_position] = mut_number #mutate parental cell

	else:#symmetric division of TAC
		state[child_position] = state[parent_position] # new TAC (on end)
		stemflag[child_position] = False
		TACAge[parent_position] += 1 #age parental TAC
		if TACAge[parent_position] == TACAgeMax:
			state = np.delete(state, parent_position)
			state = np.append(state, 0)
	return state, TACAge, stemflag, mut_pairs, mut_number


write_path = '../../figs/non_spatial/'
# cols = ['SymmDiv', 'Test', 'Value', 'iteration']
cols = ['SymmDiv', 'mutations', 'proportion', 'iteration']
df = pd.DataFrame(columns = cols)
mut_rate = 0.01
tumor_size_end = 500**2
TACAgeMax = 4

for j in range(50):
	SymDivRates = [0.2, 0.4, 0.6, 0.8, 1.0]
	for i in SymDivRates:
		state = np.zeros(tumor_size_end).astype('int')
		TACAge = np.zeros(tumor_size_end).astype('int')
		stemflag = np.zeros(tumor_size_end).astype('bool')
		state[0] = 1 # initialise
		stemflag[0] = True
		# prob_div = 0.5
		mut_number = 1
		mut_pairs = []

		while np.count_nonzero(state) < tumor_size_end:
			state, TACAge, stemflag, mut_pairs, mut_number = \
			update_non_spatial_stem_CA(state, TACAge, stemflag, mut_pairs, mut_number, i)

		# print('For case of symmetric division rate '+str(i))
		totalcells = np.count_nonzero(state)
		stemnumber = np.count_nonzero(stemflag)
		# print('Stem-proportion = ', stemnumber/(stemnumber+num_non_stems))
		# print('Total mutations = ', len(mut_pairs))

		df = df.append(pd.DataFrame([[i, len(mut_pairs), stemnumber/totalcells, j]], \
                         columns = cols), ignore_index=True)


ax = sns.boxplot(x = 'SymmDiv', y = 'mutations', data = df)
ax = sns.swarmplot(x = 'SymmDiv', y = 'mutations', data = df, color = 'k')
plt.ylabel('mutations')
# plt.show()

plt.savefig(write_path+'fixed_TAC_lowres_non_spatial_symdiv_swweep_mut'+str(mut_rate)+'.png', dpi = 250)
plt.savefig(write_path+'fixed_TAC_non_spatial_symdiv_swweep_mut'+str(mut_rate)+'.png', dpi = 500)

