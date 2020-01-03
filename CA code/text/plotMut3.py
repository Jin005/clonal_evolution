'''
Plotting CA, mutation map and mutational age.

Jacob Scott 22 Decemeber 2015

'''
import matplotlib.pyplot as plt
import numpy as np
from math import log as ln

#calculate the shannon index
def shannon(n, N):
        """ Relative abundance """
        if n is  0:
            return 0
        else:
            return (float(n)/N) * ln(float(n)/N)

#sum up the digits in a string
def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

size = 500 #size of the array
t1 = 2000
t2 = 4000
SI1 = 0 #placeholders for Shannon Index values
SI2 = 0
total_mut1 = np.zeros(size**2) #placeholders for mutation arrays
total_mut2 = np.zeros(size**2)
my_cmap = plt.cm.get_cmap('nipy_spectral')
my_cmap.set_under('w')

#timepoint 1 cells
data = open('Cells'+str(t1)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('float')
Cells1 = np.reshape(CA, (size,size))

#timepoint 1 clonal evolution
data = open('CarriedMutation'+str(t1)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + shannon(np.bincount(CA)[species],N1)

#timepoint 1 mutations
data = open('genomes'+str(t1)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('string')
for entry in range(0, size**2):	total_mut1[entry] = np.array(sum_digits(CA[entry])).astype('int')
mut_array1 = np.reshape(total_mut1, (size,size))


#timepoint 2 cells
data = open('Cells'+str(t2)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('float')
Cells2 = np.reshape(CA, (size,size))

#timepoint 2 clonal evolution
data = open('CarriedMutation'+str(t2)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM2 = np.reshape(CA, (size,size))
N2 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI2 = SI2 + shannon(np.bincount(CA)[species],N2)

#timepoint 2 mutations
data = open('genomes'+str(t2)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('string')
for entry in range(0, size**2):	total_mut2[entry] = np.array(sum_digits(CA[entry])).astype('int')
mut_array2 = np.reshape(total_mut2, (size,size))

""" PLOT """

#timepoint 1
plt.subplot(2, 3, 1)
plt.pcolor(Cells1, cmap='nipy_spectral')
plt.title('Timestep '+str(t1)+' - cells: '+str(N1))
plt.colorbar()

plt.subplot(2, 3, 2)
plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.title('Shannon Index = '+str(-SI1))
plt.colorbar()

plt.subplot(2, 3, 3)
plt.pcolor(mut_array1, cmap='nipy_spectral', vmin = 0.001)
plt.title('Total mutations')
plt.colorbar()


#timepoint 2
plt.subplot(2, 3, 4)
plt.pcolor(Cells2, cmap='nipy_spectral')
plt.title('Timestep '+str(t2)+' - cells: '+str(N2))
plt.colorbar()

plt.subplot(2, 3, 5)
plt.pcolor(CM2, cmap='nipy_spectral', vmin = 0.001)
plt.title('Shannon Index = '+str(-SI2))
plt.colorbar()

plt.subplot(2, 3, 6)
plt.pcolor(mut_array2, cmap='nipy_spectral', vmin = 0.001)
plt.title('Total mutations')
plt.colorbar()

#plt.savefig("images/3plot.png", dpi = 500)

plt.show()
