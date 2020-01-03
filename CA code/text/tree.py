'''

from ete3 import Tree
t = Tree( '((H:1,I:1):0.5, A:1, (B:1,(C:1,D:1):0.5):0.5);' )
print t
#                    /-H
#          /--------|
#         |          \-I
#         |
#---------|--A
#         |
#         |          /-B
#          \--------|
#                   |          /-C
#                    \--------|
#                              \-D

# I get D
D = t.search_nodes(name="D")[0]

# I get all nodes with distance=0.5
nodes = t.search_nodes(dist=0.5)
print len(nodes), "nodes have distance=0.5"

# We can limit the search to leaves and node names (faster method).
D = t.get_leaves_by_name(name="D")
print D

'''
import matplotlib.pyplot as plt
import numpy as np
from ete2 import Tree

data = open('output.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
ParentChild = np.array(x).astype('str')
y = len(ParentChild)/2
ParentChild1 = np.reshape(ParentChild, (y,2))
#print ParentChild1


t = Tree() # Creates an empty tree
A = t.add_child(name="A")

B = A.add_child(name="B")
C = A.add_child(name="C")
D = C.add_child(name="D")
E = A.add_child(name="E")
F = A.add_child(name="F")
G = A.add_child(name="G")
H = F.add_child(name="H") #6,8
I = A.add_child(name="I")
J = A.add_child(name="J")#10
K = D.add_child(name="K")
L = D.add_child(name="L")
M = A.add_child(name="L")#13
N = D.add_child(name="N")#4,11
O = D.add_child(name="O")#4,12
P = A.add_child(name="P")
Q = D.add_child(name="Q")#4,14
R = F.add_child(name="R")
S = D.add_child(name="S")
T = A.add_child(name="A")#1,17
U = D.add_child(name="U")#4,18
V = F.add_child(name="V")
W = A.add_child(name="W")
X = A.add_child(name="X")#1,24
Y = D.add_child(name="Y")
Z = A.add_child(name="Z")#1,26
a = A.add_child(name="AA")#1,27
BB = A.add_child(name="BB")#1,28
CC = D.add_child(name="CC")
DD = A.add_child(name="DD")
EE = A.add_child(name="EE")
FF = A.add_child(name="FF")
GG = A.add_child(name="GG")#1,33
HH = F.add_child(name="HH")#6,34
II = A.add_child(name="II")#1,35
JJ = A.add_child(name="JJ")#1,35
KK = A.add_child(name="KK")#1,35
LL = D.add_child(name="LL")#4,38
MM = A.add_child(name="MM")#1,39
NN = D.add_child(name="NN")#4,40
OO = D.add_child(name="OO")#4,41
PP = A.add_child(name="PP")#1,42
QQ = D.add_child(name="QQ")#1,39
RR = A.add_child(name="RR")#1,44
SS = F.add_child(name="SS")#6,45
TT = A.add_child(name="TT")#1,46
UU = A.add_child(name="UU")#1,47
VV = A.add_child(name="VV")#1,48
WW = A.add_child(name="WW")#1,49
XX = A.add_child(name="XX")#1,50
YY = A.add_child(name="YY")#1,51
ZZ = A.add_child(name="ZZ")#1,52
AAA = F.add_child(name="AAA")#6,53
BBB = A.add_child(name="BBB")#1,54
CCC = D.add_child(name="CCC")#4,55
DDD = A.add_child(name="DDD")#1,56
EEE = D.add_child(name="EEE")#4,57
FFF = F.add_child(name="FFF")#6,58
GGG = A.add_child(name="GGG")#1,59
HHH = A.add_child(name="HHH")#1,60
III = A.add_child(name="III")#1,61
JJJ = D.add_child(name="JJJ")#4,62
KKK = D.add_child(name="KKK")#4,63
LLL = A.add_child(name="LLL")#1,64
MMM = A.add_child(name="MMM")#1,65
OOO = D.add_child(name="OOO")#4,66
PPP = A.add_child(name="PPP")#1,67


'''

for parent in (0,y):
	c = ParentChild1[parent][1]
	p = ParentChild1[parent][0]
	print c
	c = p.add_child(name=ParentChild1[parent][1])

'''
print t.get_ascii(show_internal=True)

'''
t = Tree() # Creates an empty tree
A = t.add_child(name="A") # Adds a new child to the current tree root
                           # and returns it
B = t.add_child(name="B") # Adds a second child to the current tree
                           # root and returns it
C = A.add_child(name="C") # Adds a new child to one of the branches
D = C.add_sister(name="D") # Adds a second child to same branch as
                             # before, but using a sister as the starting
                             # point
R = A.add_child(name="R") # Adds a third child to the
                           # branch. Multifurcations are supported
# Next, I add 6 random leaves to the R branch names_library is an
# optional argument. If no names are provided, they will be generated
# randomly.
R.populate(6, names_library=["r1","r2","r3","r4","r5","r6"])
# Prints the tree topology
print t
#                     /-C
#                    |
#                    |--D
#                    |
#           /--------|                              /-r4
#          |         |                    /--------|
#          |         |          /--------|          \-r3
#          |         |         |         |
#          |         |         |          \-r5
#          |          \--------|
# ---------|                   |                    /-r6
#          |                   |          /--------|
#          |                    \--------|          \-r2
#          |                             |
#          |                              \-r1
#          |
#           \-B
# a common use of the populate method is to quickly create example
# trees from scratch. Here we create a random tree with 100 leaves.
'''


'''
t = Tree()
t.populate(100)
print t
'''