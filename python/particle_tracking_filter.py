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

# total number of particles for each species, this is only a target, 
# these will be distributed over the entire grid
Ntot_target = np.array([100000])

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


# Calculate density in center of each cell
def particle_initialization(Ntot_target,nfun):
	# 	1. calculate total density
	# 	2. loop through cells, for each cell: 
	#		2.1. calculate the average density in that cell (i just took the center density as rough estimate)
	#   	2.2. calculate the cell volume (area) dx*dy
	#		2.3. calculate how large part of the total density should go into that cell
	#		2.4. calculate how many of the total amount of particles this corresponds to
	#		2.5. calculate weight: macro-particle weight = (species density x cell hypervolume)/number of macro-particles per cell	
	n_center = np.empty([int(nx)-1, int(ny)-1])
	for i in np.arange(0,int(nx)-1):
		for j in np.arange(0,int(ny)-1):
			n_center[i,j] = nfun(x_center[i],y_center[j])

	Ntot = sum(n_center.flatten())*dx*dy # approximate total number of micro particles
	# create empty arrays
	N         = np.empty([int(nx)-1, int(ny)-1]) # N - number of micro particles in each cell
	Nfrac     = np.empty([int(nx)-1, int(ny)-1]) # N/Ntot - fraction of total number of microparticles in each cell
	M         = np.empty([int(nx)-1, int(ny)-1]) # M - number of macro particles in each cell
	x_part    = [] # particle x-coordinate
	y_part    = [] # particle y-coordinate
	weight    = [] # particle weight

	for i in np.arange(0,int(nx)-1):
		for j in np.arange(0,int(ny)-1):
			# Calculate how many macro particles should go in each cell
			cell_density = nfun(x_center[i],y_center[j])
			cell_hypervolume = dx*dy
			N[i,j] = cell_density*cell_hypervolume # density x area, number of micro particles in this cell, non-integer
			Nfrac[i,j] = N[i,j]/Ntot # fraction of micro particles that is in this cell, i.e. for a uniform density, all numbers are the same
			M[i,j] = Nfrac[i,j]*Ntot_target # number of macro particles in each cell, non-integer
			M_roundup = m.ceil(M[i,j]) # round up to nearest integer

			# Calculate weigths of particles, microparticles/macroparticles
			weight_cell = N[i,j]/M_roundup			

			# Calculate particle positions		
			xlo = x[i]
			xhi = x[i+1]
			ylo = y[j]
			yhi = y[j+1]

			new_x = np.random.uniform(xlo,xhi,size=M_roundup)
			new_y = np.random.uniform(ylo,yhi,size=M_roundup)
			new_weight = np.ones([M_roundup])*weight_cell
			x_part = np.append(x_part,new_x)
			y_part = np.append(y_part,new_y)
			weight = np.append(weight,new_weight)

	Nfinal = len(x_part)

	x_y_weight = np.empty([3,Nfinal])
	x_y_weight[0,:] = x_part
	x_y_weight[1,:] = y_part
	x_y_weight[2,:] = weight
	#print(np.shape(x_y_weight))

	#return Nfinal, x_part, y_part, weight, n_center, M

	return x_y_weight

# diagnostics
#n_part_final, x_part, y_part, weight, n_center, num_macro_particles_in_cell = particle_initialization(Ntot_target[0],n)


# run
xyw = particle_initialization(Ntot_target[0],n)
#print(np.shape(xyw))

class particles():
	def __init__(self,x,y,w):
		self.x = x
		self.y = y
		self.w = w


part = particles(xyw[0,:],xyw[1,:],xyw[2,:])

def my_filter(particles):
	return particles.x>6




if True:

	print('box_size = ' + str(box_size) + ' de')
	print('grid: ' + str(x.size) + 'x' + str(y.size) + ' cells')
	print('total number of target particles = ' + str(Ntot_target))
	#print(type(len(x)))
	#print(type(npart))
	#print('average number of particles per cell = ' + str(Ntot_target/(nx-1)/(ny-1)))
	#print('total number of particles after initiation: ' + str(Nfinal))


	# plot results


	fig, (ax1, ax2, ax2b, ax3) = plt.subplots(4 , 1)


	plt3 = ax3.scatter(part.x,part.y,s=10,c=part.w,linewidths=0,cmap='Spectral')
	fig.colorbar(plt3, ax=ax3)
	ax3.set_xlim(x[0], x[-1])
	ax3.set_ylim(y[0], y[-1])

	plt.show()

