# to be executed in terminal as >> python3 particle_weight.py
import math as m
import numpy as np
import random as rd
import matplotlib.pyplot as plt

# Here I define some python variables:
# in particular I define here the spatial & temporal resolution which do not exist anymore in the new python version
resy = 4.               # spatial resolution along y axis (i.e number of cell per de)
resx = 4.               # spatial resolution along x axis, # cells/de
rest = m.sqrt(2.*resx**2)/0.95      # time-resolution (I use dt=0.95*CFL)

# electron normalisation :

c       = 1.
wpewce  = 2.
mime    = 25.


# other parameters which I define from your notes on Harris initialisation
B0    = 1./wpewce               # amplitude of the magnetic field (we use 1 for the magnetosheat magnetic field)
n0    = 1.                      # normalization of the density
v0    = 1./(m.sqrt(mime)*wpewce)# normalization of the velocity
L     = 1.*m.sqrt(mime)        # transverse width of the B-field distribution
Lf    = 1.5*m.sqrt(mime)        # sinusoidal period of density fluctuations in inflow

#box_size = [160*m.sqrt(mime), 80*m.sqrt(mime)]
box_size = [2*m.sqrt(mime), 10*m.sqrt(mime)]
x0       = box_size[0]*0.5 #6.4 *m.sqrt(mime)        # position of the perturbation 
y0       = box_size[1]*0.5 #12.8*m.sqrt(mime)        # position of the layer

nb = 0.1*n0              # background density

# total number of particles for each species
total_num_particles = np.array([10000])

# define grid
nx, ny = list(box_size/np.array([1/resx,1/resy])+np.array([1,1]))
x = np.linspace(0,box_size[0],num=int(nx))
y = np.linspace(0,box_size[1],num=int(ny))
dx = x[1]-x[0]
dy = y[1]-y[0]
# center of grid cells, for simply calculating densities and plotting plt.pcolormesh in the center of the cell
x_center = x[0:-1] + 0.5*dx
y_center = y[0:-1] + 0.5*dy



# density profile
def n(x,y):
	return 0.2 + n0 / m.cosh((y-y0)/L)**2 
    #return nb*0.5*(1+m.tanh(((y)-2*L)/(0.5*L)))*(3+m.cos((2*np.pi/(Lf))*(y-y0)))/4

# 
# 	1. calculate total density
# 	2. loop through cells, for each cell: 
#		2.1. calculate the average density in that cell
#   	2.2. calculate the cell volume (area) dx*dy
#		2.3. calculate how large part of the total density should go into that cell
#		2.4. calculate how many of the total amount of particles this corresponds to
#		2.5. calculate weight	

# Calculate density in center of each cell
n_center = np.empty([int(nx)-1, int(ny)-1])
for i in np.arange(0,int(nx)-1):
	for j in np.arange(0,int(ny)-1):
		n_center[i,j] = n(x_center[i],y_center[j])

total_num_micro_particles = sum(n_center.flatten())*dx*dy # approximate total number of micro particles
num_particles_in_cell         = np.empty([int(nx)-1, int(ny)-1])
fraction_of_particles_in_cell = np.empty([int(nx)-1, int(ny)-1])
num_macro_particles_in_cell   = np.empty([int(nx)-1, int(ny)-1])
x_part = [] # create empty array
y_part = [] # create empty array
weight = [] # particle weight

for i in np.arange(0,int(nx)-1):
	for j in np.arange(0,int(ny)-1):
		# Calculate how many macro particles should go in each cell
		num_particles_in_cell[i,j] = n(x_center[i],y_center[j])*dx*dy # density x area, number of micro particles
		fraction_of_particles_in_cell[i,j] = num_particles_in_cell[i,j]/total_num_micro_particles
		num_macro_particles_in_cell[i,j] = fraction_of_particles_in_cell[i,j]*total_num_particles					
		
		# Calculate particle positions
		# round up to nearest integer
		# for loop probably not needed for this, just directly generate n random numbers...
		xlo = x[i]
		xhi = x[i+1]
		ylo = y[j]
		yhi = y[j+1]

		for k in np.arange(0,m.ceil(num_macro_particles_in_cell[i,j])):
			#print([i,j,k])	
			# Randomize particle position within cell
			new_x = rd.uniform(xlo, xhi)
			new_y = rd.uniform(ylo, yhi)
			x_part = np.append(x_part,new_x)
			y_part = np.append(y_part,new_y)

			# Calculate weigths of particles
			# weight = np.append(weight,...)

n_part_final = len(x_part)
#num_macro_particles_in_cell = fraction_of_particles_in_cell*total_num_particles


print('box_size = ' + str(box_size) + ' de')
print('grid: ' + str(x.size) + 'x' + str(y.size) + ' cells')
print('total number of particles = ' + str(total_num_micro_particles))
#print(type(len(x)))
#print(type(npart))
print('average number of particles per cell = ' + str(total_num_particles/(nx-1)/(ny-1)))
print('total number of particles after initiation: ' + str(n_part_final))


# plot results


fig, (ax1, ax2, ax3) = plt.subplots(3 , 1)

plt1 = ax1.pcolormesh(x,y,n_center.T,cmap='Spectral')
fig.colorbar(plt1, ax=ax1)

plt2 = ax2.pcolormesh(x,y,num_macro_particles_in_cell.T,cmap='Spectral')
fig.colorbar(plt2, ax=ax2)

plt3 = ax3.scatter(x_part,y_part,s=1)
#fig.colorbar(plt3, ax=ax3)

plt.show()

