'''
Plotting CA.

Jacob Scott 10 Decemeber 2015

'''
import matplotlib.pyplot as plt
import numpy as np
from math import log as ln
import brewer2mpl

def shannon(n, N):
        """ Relative abundance """
        if n is  0:
            return 0
        else:
            return (float(n)/N) * ln(float(n)/N)

#bmap = brewer2mpl.get_map('Paired', 'Qualitative', 5, reverse=True)

size = 1000 #size of the array
SI1 = 0
SI2 = 0
my_cmap = plt.cm.get_cmap('nipy_spectral')
my_cmap.set_under('w')

#newcolormap = [.RdYlBu .q0-3{fill:rgb(252,141,89)} .RdYlBu .q1-3{fill:rgb(255,255,191)} .RdYlBu .q2-3{fill:rgb(145,191,219)}]

#timepoint 1 cells
data = open('cells1000').read().replace(',','\n').replace('\n\n','\n').replace('\n','').replace('\n','').replace('  ',' ')
data = data[:-1]
x = data.split()
CA = np.array(x).astype('float')
CA1 = np.reshape(CA, (size,size))

#timepoint 1 clonal evolution
data = open('carriedmutation1000').read().replace(',','\n').replace('\n\n','\n').replace('\n','').replace('\n','').replace('  ',' ')
data = data[:-1]
x = data.split()
CA = np.array(x).astype('int')
CA2 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + shannon(np.bincount(CA)[species],N1)

#timepoint 2 cells
data = open('cells2000').read().replace(',','\n').replace('\n\n','\n').replace('\n','').replace('\n','').replace('  ',' ')
data = data[:-1]
x = data.split()
CA = np.array(x).astype('float')
CA3 = np.reshape(CA, (size,size))

#timepoint 2 clonal evolution
data = open('carriedmutation2000').read().replace(',','\n').replace('\n\n','\n').replace('\n','').replace('\n','').replace('  ',' ')
data = data[:-1]
x = data.split()
CA = np.array(x).astype('int')
CA4 = np.reshape(CA, (size,size))
N2 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI2 = SI2 + shannon(np.bincount(CA)[species],N2)

#print CA
plt.subplot(2, 2, 1)
plt.pcolor(CA1, cmap='brg')
plt.title('Timestep 1000 - cells: '+str(N1))
plt.colorbar()

plt.subplot(2, 2, 2)
plt.pcolor(CA2, cmap='nipy_spectral', vmin = 0.001)
plt.title('Shannon Index = '+str(-SI1))
plt.colorbar()

plt.subplot(2, 2, 3)
plt.pcolor(CA3, cmap='brg')
plt.title('Timestep 2000 - cells: '+str(N2))
plt.colorbar()

ax = plt.subplot(2, 2, 4)
plt.pcolor(CA4, cmap='nipy_spectral', vmin = 0.001)
plt.title('Shannon Index = '+str(-SI2))
plt.colorbar()

plt.show()
