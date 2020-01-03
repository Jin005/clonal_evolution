'''
non-spatial CA to serve as null model to study effects of spatial organization on phylogenies

september 21st, 2016

jacob g. scott
'''
import dendropy
from dendropy import Tree as DTree
from dendropy.calculate import treemeasure
import numpy as np
import matplotlib.pyplot as plt
from random import random
from ete3 import Tree
import clonal_evolution_functions as cef
import pandas as pd
import seaborn as sns
from pylab import rcParams
import CSC_functions as csc

write_path = '../../figs/non_spatial/'
cols = ['SymmDiv', 'Test', 'Value', 'iteration']
cols_colonies = ['SymmDiv', 'Test', 'Value', 'iteration']
# cols = ['SymmDiv', 'mutations', 'proportion', 'iteration']
df = pd.DataFrame(columns = cols)
df_colonies = pd.DataFrame(columns = cols_colonies)
mut_rate = 0.01
tumor_size_end = 500**2
TACAgeMax = 4

for j in range(1):
	SymDivRates = [0.2, 0.4, 0.6, 0.8, 1.0]
	for i in SymDivRates:
		state = np.zeros(tumor_size_end).astype('int')
		TACAge = np.zeros(tumor_size_end).astype('int')
		stemflag = np.zeros(tumor_size_end).astype('bool')
		state[0] = 1 # initialise
		stemflag[0] = True
		mut_number = 1
		mut_pairs = []

		while np.count_nonzero(state) < tumor_size_end:
			state, TACAge, stemflag, mut_pairs, mut_number = \
			csc.update_non_spatial_stem_CA(state, TACAge, TACAgeMax, stemflag, mut_pairs, mut_number, mut_rate, i)

###########


		ete_Tree = cef.make_tree_from_list(mut_pairs)
		
		# print(ete_Tree)

		treePass = ete_Tree.write(format = 1)
		DendroTree = DTree.get(data = treePass, schema = 'newick')


		''' Tree is now in dendropy format and ready to measure '''

		df = df.append(pd.DataFrame([[i, 'B1', treemeasure.B1(DendroTree), j]], \
		                         columns = cols), ignore_index=True)

		df = df.append(pd.DataFrame([[i, 'Nbar', treemeasure.N_bar(DendroTree), j]], \
		                         columns = cols), ignore_index=True)

		df = df.append(pd.DataFrame([[i, 'yule', treemeasure.sackin_index(DendroTree, normalize = 'yule'), j]], \
		                         columns = cols), ignore_index=True)

		df = df.append(pd.DataFrame([[i, 'PDA', treemeasure.sackin_index(DendroTree, normalize = 'pda'), j]], \
		                         columns = cols), ignore_index=True)

		df = df.append(pd.DataFrame([[i, 'Sackin', treemeasure.sackin_index(DendroTree, normalize = None), j]], \
		                         columns = cols), ignore_index=True)

		# print(sizes[i],j)

# df.to_csv(write_path+'multi_treeMetric_non_spatial_highsize_'+str(mut_rate)+'.csv')

rcParams['figure.figsize'] = 12, 12

plt.subplot(221)
df_yule = df[df['Test'] == 'yule']
ax = sns.boxplot(x = 'SymmDiv', y = 'Value', data = df_yule)
ax = sns.swarmplot(x = 'SymmDiv', y = 'Value', data = df_yule, color = 'k')
plt.ylabel('Yule')

plt.subplot(222)
df_Nbar = df[df['Test'] == 'Nbar']
ax = sns.boxplot(x = 'SymmDiv', y = 'Value', data = df_Nbar)
ax = sns.swarmplot(x = 'SymmDiv', y = 'Value', data = df_Nbar, color = 'k')
plt.ylabel('Nbar')

plt.subplot(223)
df_sackin = df[df['Test'] == 'Sackin']
ax = sns.boxplot(x = 'SymmDiv', y = 'Value', data = df_sackin)
ax = sns.swarmplot(x = 'SymmDiv', y = 'Value', data = df_sackin, color = 'k')
plt.ylabel('Sackin')

plt.subplot(224)
df_B1 = df[df['Test'] == 'B1']
ax = sns.boxplot(x = 'SymmDiv', y = 'Value', data = df_B1)
ax = sns.swarmplot(x = 'SymmDiv', y = 'Value', data = df_B1, color = 'k')
plt.ylabel('B1')
# print(t)
plt.savefig('TIMETESTnon_spatial_trees_low_res'+str(mut_rate)+'.png', dpi = 250)
plt.savefig('TIMETESTnon_spatial_trees'+str(mut_rate)+'.png', dpi = 500)
# plt.show()


		# totalcells = np.count_nonzero(state)
		# stemnumber = np.count_nonzero(stemflag)

		# df = df.append(pd.DataFrame([[i, len(mut_pairs), stemnumber/totalcells, j]], \
  #                        columns = cols), ignore_index=True)


# ax = sns.boxplot(x = 'SymmDiv', y = 'mutations', data = df)
# ax = sns.swarmplot(x = 'SymmDiv', y = 'mutations', data = df, color = 'k')
# plt.ylabel('mutations')
# plt.show()

# plt.savefig(write_path+'full_lowres_non_spatial_symdiv_swweep_mut'+str(mut_rate)+'.png', dpi = 250)
# plt.savefig(write_path+'full_non_spatial_symdiv_swweep_mut'+str(mut_rate)+'.png', dpi = 500)

# plt.savefig(write_path+'test_full'+str(mut_rate)+'.png', dpi = 250)



