''' 
Build multiple, temporally distinct phylogenies from (parent,child) lists

Jacob G Scott, 31 May 2016 
'''


import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace
from collections import Counter
from math import log as ln
import random as random
import sys


# read_path = str(sys.argv[1]) 
read_path = '../patchiness_experiments/dot4/multi_sim_tree_outputs/'
# write_path = str(sys.argv[2])
write_path = '../patchiness_experiments/dot4/figs/time_loops/'
# filename = str(sys.argv[3])
filename = 'multi_sim_output_dot4_1'

def sort_pairs(pair):
    # Extract integer after "r".
    return int(pair[0][1:])

# def pause():
#     programPause = input("Press the <ENTER> key to continue...")

# #life history information
# data = open(read_path+filename).read().replace(',',' ').replace('\n',' ')
# x = data.split()
# x = x[3:] #cuts off first 3 entries which are simply the simulation parameters
# ParentChild = np.array(x).astype(str)
# y = len(ParentChild)/5
# ParentChild1 = np.reshape(ParentChild, (y,5))
# firsttwo = ParentChild1[:,0:2] #chops off columns 3,4,5 which are stem, non-stem, and time


#routine that makes N version of firsttwo each starting at the beginning but getting longer by 1/N
N = 20

# print(firsttwo)
for i in range(1,N+1):
#life history information
    data = open(read_path+filename).read().replace(',',' ').replace('\n',' ')
    x = data.split()
    x = x[3:] #cuts off first 3 entries which are simply the simulation parameters
    ParentChild = np.array(x).astype(str)
    y = len(ParentChild)/5
    ParentChild1 = np.reshape(ParentChild, (y,5))
    firsttwo = ParentChild1[:,0:2] #chops off columns 3,4,5 which are stem, non-stem, and time
    
    history = firsttwo[0:i*len(firsttwo)/(N)]
    # print('firsttwo: ',firsttwo[0:i*len(firsttwo)/(N)])
    # print('history: ',history)

    parents = []
    children = []

    prune_list = []
    lookup = {}

    for row in range(0, len(history)): 
    	for column in range(0,2): 
    		history[row,column] = 'r'+history[row,column]
    # print('rs: ',history)

    t = Tree() # Creates an empty tree

    r1 = t.add_child(name="r1")
    lookup = {"r1": r1}
    prune_list = ['r1']

    for pair in sorted(history, key=sort_pairs):
        parentname = pair[0]
        childname = pair[1]
        if childname not in lookup:
            if parentname in lookup:
                newchild = lookup[parentname].add_child(name = childname)
                lookup.update({childname: newchild})
                if parentname not in parents:
                    prune_list.append(lookup[parentname])
                parents.append(parentname) #make list of unique terminal nodes (no children of children)
                children.append(newchild)
            else:
                raise RuntimeError('Must not happen.')
    # print('history: ',history)
    # pause()

    # prune_count = Counter(children) #counter than contains the number of children that each terminal node has]


    t.write(outfile = write_path+'unpruned_history_'+str(i)+'_'+filename+'.nw', format = 100)

    t.prune(prune_list)

    # print (t.get_ascii(show_internal=True))
    t.write(outfile = write_path+'pruned_history_'+str(i)+'_'+filename+'.nw', format = 100)
    # t.show(tree_style=ts)
    