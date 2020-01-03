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
# import seaborn as sns
from pylab import rcParams


"""discrete population dynamics, stem"""
def update_non_spatial_stem_CA(state, mut_pairs, mut_number, SymDiv):
	parent = state[np.random.randint(np.count_nonzero(state))].astype(int) #choose cell to divide
	num_non_stems = np.bincount(state)[tumor_size_end+1] # count up non-stems
	if parent < tumor_size_end+1: #if chosen cell is a stem cell
		if rand() < SymDiv: #is this symmetric division?
			num_muts = np.random.poisson(mut_rate, 1).astype(int) #get some mutations
			if num_muts == 0:
				state[np.count_nonzero(state)] = parent # new daughter (Tacked on end), is unmutated
			else:
				parentname = parent.astype(str)
				child = (mut_number + num_muts).astype(int)
				childname = child[0].astype(str)
				mut_pairs.append((parentname, childname))
				mut_number = mut_number+num_muts
				state[np.count_nonzero(state)+num_non_stems] = mut_number #add mutated cell at the end
		else: 	state[np.count_nonzero(state)+num_non_stems] = tumor_size_end+1 # new TAC (on end)
	else: 	state[np.count_nonzero(state)+num_non_stems] = tumor_size_end+1 # new TAC (on end)
	return state, mut_pairs, mut_number


write_path = '../figs/non_spatial/'
cols = ['SymmDiv', 'Test', 'Value', 'iteration']
df = pd.DataFrame(columns = cols)
mut_rate = 0.05
tumor_size_end = 5**2

for j in range(20):
	SymDiv = [0.2, 1.0]
	for i in range(len(SymDiv)):
		state = np.zeros(tumor_size_end)
		state[0] = 1 # initialise
		time = 0
		prob_div = 1
		mut_number = 1
		mut_pairs = []

		# print(state)
		# print(np.count_nonzero(state))

		while np.count_nonzero(state) < tumor_size_end:
			state, mut_pairs, mut_number = update_non_spatial_stem_CA(state, mut_pairs, mut_number, i)
			# print(state)
		# print('All mutational pairs', mut_pairs)
		# print('Final state vector', state)
		# print('Tumor size = ', np.count_nonzero(state))


############


		# ete_Tree = cef.make_tree_from_list(mut_pairs)
		
		# print(ete_Tree)

# 		treePass = ete_Tree.write(format = 1)
# 		DendroTree = DTree.get(data = treePass, schema = 'newick')

# 		''' Tree is now in dendropy format and ready to measure '''

# 		df = df.append(pd.DataFrame([[sizes[i], 'B1', treemeasure.B1(DendroTree), j]], \
# 		                         columns = cols), ignore_index=True)

# 		df = df.append(pd.DataFrame([[sizes[i], 'Nbar', treemeasure.N_bar(DendroTree), j]], \
# 		                         columns = cols), ignore_index=True)

# 		df = df.append(pd.DataFrame([[sizes[i], 'yule', treemeasure.sackin_index(DendroTree, normalize = 'yule'), j]], \
# 		                         columns = cols), ignore_index=True)

# 		df = df.append(pd.DataFrame([[sizes[i], 'PDA', treemeasure.sackin_index(DendroTree, normalize = 'pda'), j]], \
# 		                         columns = cols), ignore_index=True)

# 		df = df.append(pd.DataFrame([[sizes[i], 'Sackin', treemeasure.sackin_index(DendroTree, normalize = None), j]], \
# 		                         columns = cols), ignore_index=True)

# 		print(sizes[i],j)

# df.to_csv(write_path+'multi_treeMetric_non_spatial_highsize_mut05.csv')

# rcParams['figure.figsize'] = 15, 15

# plt.subplot(221)
# df_yule = df[df['Test'] == 'yule']
# ax = sns.boxplot(x = 'Size', y = 'Value', data = df_yule)
# ax = sns.swarmplot(x = 'Size', y = 'Value', data = df_yule, color = 'k')
# plt.ylabel('Yule')

# plt.subplot(222)
# df_Nbar = df[df['Test'] == 'Nbar']
# ax = sns.boxplot(x = 'Size', y = 'Value', data = df_Nbar)
# ax = sns.swarmplot(x = 'Size', y = 'Value', data = df_Nbar, color = 'k')
# plt.ylabel('Nbar')

# plt.subplot(223)
# df_sackin = df[df['Test'] == 'Sackin']
# ax = sns.boxplot(x = 'Size', y = 'Value', data = df_sackin)
# ax = sns.swarmplot(x = 'Size', y = 'Value', data = df_sackin, color = 'k')
# plt.ylabel('Sackin')

# plt.subplot(224)
# df_B1 = df[df['Test'] == 'B1']
# ax = sns.boxplot(x = 'Size', y = 'Value', data = df_B1)
# ax = sns.swarmplot(x = 'Size', y = 'Value', data = df_B1, color = 'k')
# plt.ylabel('B1')
# # print(t)
# plt.savefig(write_path+'4plot_non_spatial_stem_sweep_lowres_mut'+str(mut_rate)+'.png', dpi = 250)
# plt.savefig(write_path+'4plot_non_spatial_stem_sweep_mut'+str(mut_rate)+'.png', dpi = 500)
# # plt.show()

