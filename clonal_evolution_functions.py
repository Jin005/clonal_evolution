''' Functions to be used for analysis of clonal evolution CA 

Jacob G Scott Jan 2016

'''
from math import log as ln
from math import ceil as ceil
from math import sqrt as sqrt
import random
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import scipy.stats as stats
from math import ceil as ceil
from ete3 import Tree


def pause():
    programPause = input("Press the <ENTER> key to continue...")

def pdf(name='temp', show=True):
    """figsave name -> pylab.savefig('name.pdf', transparent = True)"""
    f = plt.gcf()
    for ax in f.axes:
        ax.de_clip(justTicks=True)
    name += '.pdf'
        
    plt.savefig(name, format='pdf', trransparent=True)
    if show:
        Popen(['xdg-open' if sys.platform == 'linux' else 'open', name])

def sort_pairs(pair):
    # Extract integer after "r".
    return int(pair[0][1:])

def make_tree_from_list(mut_pairs):
	parents = []
	children = []
	pairs_of_mutations = []
	for item in mut_pairs:
		a = 'r'+item[0]
		b = 'r'+item[1]
		pairs_of_mutations.append((a,b))
	t = Tree() # Creates an empty tree
	r1 = t.add_child(name="r1")
	lookup = {"r1": r1}
	for pair in sorted(pairs_of_mutations, key=sort_pairs):
		parentname = pair[0]
		childname = pair[1]
		if childname not in lookup:
		    if parentname in lookup:
		        newchild = lookup[parentname].add_child(name = childname)
		        lookup.update({childname: newchild})
		        parents.append(parentname) #make list of unique terminal nodes (no children of children)
		        children.append(newchild)
		    else:
		        raise RuntimeError('Must not happen.')
	return t		        

''' calculate the shannon index '''
def shannon(n, N):
    if n == 0:
        return 0
    else:
        return (float(n)/N) * ln(float(n)/N)

''' sum up the digits in a string '''
def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

''' Function to count mutations in entire array ''' 
def total_mutation_map(CM1, size, family_dict):
	mutation_number = np.zeros((size,size)).astype(int)
	for (row, col), cell in np.ndenumerate(CM1):
		if cell == 0: continue
		mutation_number[row][col] += 1
		if cell == 1: continue
		mutation_number[row][col] += 1
		parent = family_dict[cell]
		while parent > 1:
			mutation_number[row][col] += 1
			parent = family_dict[parent]
	return mutation_number

'''find neighbors of prescribed size for biopsy'''
def neighborhood(biopsy_site, r):
	neighboring_cells = []
	x = biopsy_site[0]
	y = biopsy_site[1]
	for i in range(0,r+1):
		for j in range(0,r+1):
			if distance.euclidean((x+i,y+j),biopsy_site) <= r and (x+i,y+j) not in neighboring_cells:
				neighboring_cells.append((x+i, y+j))
			if distance.euclidean((x-i,y+j),biopsy_site) <= r and (x-i,y+j) not in neighboring_cells:
				neighboring_cells.append((x-i, y+j))
			if distance.euclidean((x-i,y-j),biopsy_site) <= r and (x-i,y-j) not in neighboring_cells:
				neighboring_cells.append((x-i, y-j))
			if distance.euclidean((x+i,y-j),biopsy_site) <= r and (x+i,y-j) not in neighboring_cells:
				neighboring_cells.append((x+i, y-j))
	return neighboring_cells #returns ordered pairs

''' Take biopsies '''
def do_biopsies(size, biopsy_num, r, CM1, biopsy_sites):
	area = 4*r**2
	biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
	cell_count_inBx = np.zeros(biopsy_num)
	SIBx = []
	for row in range(0, size):
		for column in range(0, size):
			a = (row,column)
			for bx in range(0, biopsy_num):
				punch = biopsy_sites[bx]
				if distance.euclidean(a,punch) <= r:
					biopsy_Mutlist[bx][cell_count_inBx[bx]] = CM1[column][row]
					cell_count_inBx[bx] += 1
	for bx in range(0, biopsy_num):
		SIBx_temp = 0
		biopsy_Mutlist_temp = (biopsy_Mutlist[bx])[0:cell_count_inBx[bx]]
		for x in range (0, np.amax(biopsy_Mutlist_temp)):
			SIBx_temp += shannon(np.bincount(biopsy_Mutlist_temp)[x],cell_count_inBx[bx])
		SIBx_temp = float("{0:.3f}".format(SIBx_temp))
		SIBx.append(-SIBx_temp)
	return SIBx

''' Function to derive genome and count mutations in provided list of cells ''' 
def derive_genome_biopsy(biopsy_list, family_dict, CM1):
	derived_genomes_inBx = list(map(str, np.zeros(len(biopsy_list))))
	for biopsy in range(0,len(biopsy_list)):
		if biopsy_list[biopsy] == 0:
			bitstring = (np.max(CM1))*'0'
			derived_genomes_inBx[biopsy] = ''.join(bitstring)
			continue
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if biopsy_list[biopsy] == 1:
			derived_genomes_inBx[biopsy] = ''.join(bitstring)
			continue 
		else:
			temp_parent = family_dict[biopsy_list[biopsy]]
			bitstring[biopsy_list[biopsy]-1] = '1'
			while temp_parent > 1:
				temp_parent = family_dict[temp_parent]
				bitstring[temp_parent-1] = '1'
				if temp_parent == 1: break			
			derived_genomes_inBx[biopsy] = ''.join(bitstring)
	return derived_genomes_inBx

''' Function to count mutations in entire array ''' 
def count_mutations(CM1, family_dict):
	mutation_number = np.zeros((size,size)).astype(int)
	for (row, col), cell in np.ndenumerate(CM1):
		if cell == 0: continue
		temp_parent = 2
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if cell == 1:
			mutation_number[row][col] = sum_digits(bitstring)
			continue 
		else:
			while temp_parent > 1:
				temp_parent = family_dict[cell]
				bitstring[cell-1] = '1'
				if temp_parent == 1: break
				cell = family_dict[cell]
			mutation_number[row][col] = sum_digits(bitstring)
	return mutation_number


''' Function to count mutations and derive genome in entire array ''' 
def count_derive_mutations(CM1, family_dict):
	mutation_number = np.zeros((size,size)).astype(int)
	derived_genomes = np.zeros(len(biopsy_list)).astype(str)
	for (row, col), cell in np.ndenumerate(CM1):
		if cell == 0: continue
		temp_parent = 2
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if cell == 1:
			derived_genomes[row][col] = ''.join(bitstring)
			mutation_number[row][col] = sum_digits(bitstring)
			continue 
		else:
			while temp_parent > 1:
				temp_parent = family_dict[cell]
				bitstring[cell-1] = '1'
				if temp_parent == 1: break
				cell = family_dict[cell]
			derived_genomes[row][col] = ''.join(bitstring)
			mutation_number[row][col] = sum_digits(bitstring)
	return mutation_number,derived_genomes


''' Identify centers of areas to biopsy and return them, r currently redundant '''
def gather_biopsies(biopsy_num, r, size):
	biopsy_sites = []
	point1 = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #pick a random position at least r from the edge
	biopsy_sites.append(point1)
	while len(biopsy_sites) < biopsy_num:
		newpoint = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #not including over the edge
		# newpoint = [random.randint(0,size),random.randint(0,size)] #overlap OK, over edge ok
		biopsy_sites.append(newpoint)
	return biopsy_sites

''' Gather biopsies from square, diamond or both (middle 1/4 and at edge 1/8) '''
def square_biopsies(size, mid_square, mid_diamond, mid_circ, edge_square, edge_diamond, edge_circ):
	biopsy_sites = []
	if mid_square == True: # at 1/4 and 3/4
		biopsy_sites.append((ceil(size/4), ceil(size/4)))
		biopsy_sites.append((ceil(3*size/4), ceil(size/4)))
		biopsy_sites.append((ceil(size/4), ceil(3*size/4)))
		biopsy_sites.append((ceil(3*size/4), ceil(3*size/4)))
	if edge_square == True: # at 1/8 and 7/8
		biopsy_sites.append((ceil(size/8), ceil(size/8)))
		biopsy_sites.append((ceil(7*size/8), ceil(size/8)))
		biopsy_sites.append((ceil(size/8), ceil(7*size/8)))
		biopsy_sites.append((ceil(7*size/8), ceil(7*size/8)))
	if mid_circ == True: # at 1/4 and 3/4
		biopsy_sites.append((ceil(size/2 + size/(4*sqrt(2))), ceil(size/2 + size/(4*sqrt(2)))))
		biopsy_sites.append((ceil(size/2 - size/(4*sqrt(2))), ceil(size/2 + size/(4*sqrt(2)))))
		biopsy_sites.append((ceil(size/2 + size/(4*sqrt(2))), ceil(size/2 - size/(4*sqrt(2)))))
		biopsy_sites.append((ceil(size/2 - size/(4*sqrt(2))), ceil(size/2 - size/(4*sqrt(2)))))
	if edge_circ == True: # at 1/8 and 7/8
		biopsy_sites.append((ceil(size/2 + 3*size/(8*sqrt(2))), ceil(size/2 + 3*size/(8*sqrt(2)))))
		biopsy_sites.append((ceil(size/2 + 3*size/(8*sqrt(2))), ceil(size/2 - 3*size/(8*sqrt(2)))))
		biopsy_sites.append((ceil(size/2 - 3*size/(8*sqrt(2))), ceil(size/2 + 3*size/(8*sqrt(2)))))
		biopsy_sites.append((ceil(size/2 - 3*size/(8*sqrt(2))), ceil(size/2 - 3*size/(8*sqrt(2)))))
	if mid_diamond == True:
		biopsy_sites.append((ceil(size/4), ceil(size/2)))
		biopsy_sites.append((ceil(3*size/4), ceil(size/2)))
		biopsy_sites.append((ceil(size/2), ceil(3*size/4)))
		biopsy_sites.append((ceil(size/2), ceil(size/4)))
	if edge_diamond == True:
		biopsy_sites.append((ceil(size/8), ceil(size/2)))
		biopsy_sites.append((ceil(7*size/8), ceil(size/2)))
		biopsy_sites.append((ceil(size/2), ceil(7*size/8)))
		biopsy_sites.append((ceil(size/2), ceil(size/8)))	
	return biopsy_sites

''' Identify centers of areas to biopsy and return them, r currently redundant '''
def gather_biopsies_noOverlap(biopsy_num, r, size):
	biopsy_sites = []
	point1 = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #pick a random position at least r from the edge
	biopsy_sites.append(point1)
	while len(biopsy_sites) < biopsy_num:
		distance_list = []
		newpoint = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #not including over the edge
		for site in biopsy_sites:
			distance_list.append(distance.euclidean(site, newpoint))
		if min(distance_list) >= 2*r: biopsy_sites.append(newpoint)	
	return biopsy_sites

''' Gathers a given number of biopsies and reports the distance from one another '''
def gather_spaced_biopsies(r, size, equidistant = False, biopsy_num = 2):
	biopsy_sites = []
	point1 = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #pick a random position at least r from the edge
	while distance.euclidean(point1,(size/2, size/2)) < r:
			point1 = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #pick a random position at least r from the edge
	biopsy_sites.append(point1)
	while len(biopsy_sites) < biopsy_num:
		newpoint = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #not including over the edge
		if equidistant == True: 
			if (distance.euclidean(biopsy_sites[0], (size/2,size/2)) - r/2) < distance.euclidean(newpoint, (size/2, size/2)) < (distance.euclidean(biopsy_sites[0], (size/2,size/2)) + r/2):
				biopsy_sites.append(newpoint)
				ptptdist = distance.euclidean(biopsy_sites[0],newpoint)
		elif distance.euclidean(newpoint, biopsy_sites[0]) > 2*r:
			ptptdist = distance.euclidean(biopsy_sites[0],newpoint)
			biopsy_sites.append(newpoint)
	ptcenterdist = distance.euclidean(newpoint, ((size/2),(size/2))) + distance.euclidean(point1, ((size/2),(size/2)))
	return biopsy_sites, ptptdist, ptcenterdist

''' Gathers pairs of biopsies the same distance from the center, but reflected through the origin'''
def gather_reflected_biopsies(size, r, center_distance):
	biopsy_sites = []
	while len(biopsy_sites) < 1:
		point1 = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #pick a random position at least r from the edge
		if (center_distance - r/2) < distance.euclidean(point1, (size/2,size/2)) < (center_distance + r/2):
			biopsy_sites.append(point1)
	x2 = biopsy_sites[0][0] + 2*(size/2 - biopsy_sites[0][0])
	y2 = biopsy_sites[0][1] + 2*(size/2 - biopsy_sites[0][1])
	biopsy_sites.append((x2,y2))
	return biopsy_sites

''' Take biopsies and return a list of the mutations present and number of cells '''
def return_biopsied_mutations(size, biopsy_num, r, CM1, biopsy_sites):
	area = 4*r**2
	biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
	cell_count_inBx = np.zeros(biopsy_num)
	for (row, col), cell in np.ndenumerate(CM1):
		a = (row,col)
		for bx in range(0, biopsy_num):
			punch = biopsy_sites[bx]
			if distance.euclidean(a,punch) <= r:
				biopsy_Mutlist[bx][cell_count_inBx[bx]] = cell
				cell_count_inBx[bx] += 1
	return biopsy_Mutlist, cell_count_inBx

''' Take biopsies and return a list of the mutations present and number of cells '''
def return_bx_muts_from_sites(CM1, neighboring_cells):
	biopsy_Mutlist = np.zeros(len(neighboring_cells)).astype('int')
	for i in range (0,len(neighboring_cells)):
		pair = neighboring_cells[i]	
		biopsy_Mutlist[i] = CM1[pair[1]][pair[0]]
	return biopsy_Mutlist

''' Gather biopsies from CA, aggregate all cells into one list and calculate SI ''' 
def agg_bx_fast(biopsy_num, CM1, r, biopsy_sites):
	SIBx_agg = 0
	SIBx = []
	biopsied_mutations = []
	aggregate_biopsy = np.array([]).astype('int')
	for i in range(0,biopsy_num):
		SIBx_temp = 0
		neighborhood_temp = neighborhood(biopsy_sites[i],r)
		biopsied_mutations = return_bx_muts_from_sites(CM1, neighborhood_temp)
		for x in range (0,np.amax(biopsied_mutations)):
			SIBx_temp += shannon(np.bincount(biopsied_mutations)[x],len(neighborhood_temp))
		SIBx_temp = float("{0:.3f}".format(SIBx_temp))
		SIBx.append(-SIBx_temp)
		aggregate_biopsy = np.append(aggregate_biopsy, biopsied_mutations)
	for x in range (0, np.amax(aggregate_biopsy)):
		SIBx_agg += shannon(np.bincount(aggregate_biopsy)[x],biopsy_num*len(neighborhood_temp))
	SIBx_agg = float("{0:.3f}".format(-SIBx_agg))
	return SIBx, SIBx_agg



''' OLD FUNCTIONS '''
''' Gather biopsies from CA, aggregate all cells into one list and calculate SI ''' 
# def do_biopsies_aggregate(size, biopsy_num, r, CM1, biopsy_sites):
# 	area = 4*r**2
# 	biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
# 	aggregate_biopsy = np.array([]).astype('int')
# 	cell_count_inBx = np.zeros(biopsy_num)
# 	SIBx_agg = 0
# 	for row in range(0, size):
# 		for column in range(0, size):
# 			a = (row,column)
# 			for bx in range(0, biopsy_num):
# 				punch = biopsy_sites[bx]
# 				if distance.euclidean(a,punch) <= r:
# 					biopsy_Mutlist[bx][cell_count_inBx[bx]] = CM1[column][row]
# 					cell_count_inBx[bx] += 1 
# 					aggregate_biopsy = np.append(aggregate_biopsy,CM1[column][row])
# 	for bx in range(0, biopsy_num):
# 		SIBx_temp = 0
# 		biopsy_Mutlist_temp = (biopsy_Mutlist[bx])[0:cell_count_inBx[bx]]
# 		for x in range (0, np.amax(biopsy_Mutlist_temp)):
# 			SIBx_temp += shannon(np.bincount(biopsy_Mutlist_temp)[x],cell_count_inBx[bx])
# 		SIBx_temp = float("{0:.3f}".format(SIBx_temp))
# 		SIBx.append(-SIBx_temp)
# 	for x in range (0, np.amax(aggregate_biopsy)):
# 		SIBx_agg += shannon(np.bincount(aggregate_biopsy)[x],np.sum(cell_count_inBx))
# 	SIBx_agg = float("{0:.3f}".format(SIBx_agg))
# 	return SIBx, -SIBx_agg

