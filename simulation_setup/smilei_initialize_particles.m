% Smilei accepts hdf5 files as particle input.
% This is useful for dealing with an a skew distribution of number of macro
% particles per microparticles. I.e. the Harris sheet has zero particles
% away from the current sheet. 
% ----------------------------------------------------------------------
% https://smileipic.github.io/Smilei/particle_initialization.html
% The position_initialization may be a path to an HDF5 file containing the 
% appropriate data structure. The path may point to a file, such as 
% "some_folder/some_data.h5", but it may also contain the path to a group 
% inside the file, such as "some_folder/some_data.h5/group1/group2".
% 
% The HDF5 location must contain the following datasets, all 1-dimensional of equal size:
% 
%     position/x, list of x coordinates
% 
%     position/y, list of y coordinates
% 
%     position/z, list of z coordinates
% 
%     weight, list of statistical weights
% 
% The momentum_initialization works the same way. It must be an HDF5 
% location containing the following datasets, all 1-dimensional of equal 
% size:
% 
%     momentum/x, list of px
% 
%     momentum/y, list of py
% 
%     momentum/z, list of pz
% 
% Note
% 
% This file structure is identical to that obtained from the TrackParticles
% diagnostics, meaning that you can directly pass the output of a previous 
% simulation, for instance 
% "path/to/results/TrackParticlesDisordered_myspecies.h5/data/0000003000/particles/myspecies".

% Define agrid into which you divide grid. The grid does not have to be the
% same as in the simulation, since it's only used to discretize n. But xmax
% and xmin needs to be the same as in simulation.
xmin = 0; xmax = 128/2;
ymin = 0; ymax = 64/2;
nx = 100;6401; % edges points of grid in z
ny = 50;1281; % edges points of grid in y
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
y0 = mean(y);
dx = x(2) - x(1);
dy = y(2) - y(1);
cell_volume = dx*dy;
x_center = x(2:end)-0.5*dx;
y_center = y(2:end)-0.5*dy;

% This is the shape of your density distribution
n0 = 1;
L = 2;
%nfun = @(x,y) x*0 + 1;
nfun = @(x,y) 0.001 + n0./cosh((y-y0)/L).^2; 

% Target number of macro_particles int he simulation (for given species)
M_tot_target = 1e5;

% Discretize density distribution onto centers of grid.
[X,Y] = ndgrid(x_center,y_center);
N = nfun(X,Y);

% Total number of microparticles in simulation (for given species)
N_tot = sum(N(:)); 

% Fraction of particles in each cell
N_frac = N/N_tot;

% Target number of macro particles in each cell
M_target = M_tot_target*N_frac;

% Final number of macro particles in each cell, round up M_target
M = ceil(M_target);

% Final particle weight, microparticles/macroparticles
W = N./M;

% Final total number of all macro particles
M_final = sum(M(:));

% Distribute the assigned number of particles throughout the cell
% It would be best to do it in some matrix form.
[X_left_edges,Y_left_edges] = ndgrid(x(1:end-1),y(1:end-1));

% Initialize matrix for particle positions
part_x = zeros(M_final,1);
part_y = zeros(M_final,1);
part_w = zeros(M_final,1);
% Fill matrix with random positions within each cell
part_count_stop = 0;
for icell = 1:numel(X_left_edges)
  icell;
  part_count_start = part_count_stop  + 1;
  part_count_stop  = part_count_start - 1 + M(icell);
  
  part_x(part_count_start:part_count_stop) = repmat(X_left_edges(icell),M(icell),1) + dx*rand(M(icell),1);
  part_y(part_count_start:part_count_stop) = repmat(Y_left_edges(icell),M(icell),1) + dy*rand(M(icell),1);
  part_w(part_count_start:part_count_stop) = W(icell);
end

% Assign thermal velocities for each particle

vt = 0.1;




% Plot results
nrows = 4;
ncols = 1;

h = setup_subplots(nrows,ncols);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,X,Y,N)
shading(hca,'flat');
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'N';


hca = h(isub); isub = isub + 1;
pcolor(hca,X,Y,N_frac)
shading(hca,'flat');
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Nfrac';


hca = h(isub); isub = isub + 1;
pcolor(hca,X,Y,M_target)
shading(hca,'flat');
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Mtarget';

hca = h(isub); isub = isub + 1;
plot(hca,part_x,part_y,'.')
%shading(hca,'flat');
%hcb = colorbar('peer',hca);
%hcb.YLabel.String = 'Mtarget';
% PYTHON CODE
% ----------------------------------------------------------------------
% def particle_initialization(Mtot_target,nfun):    
%     #   1. loop through cells, calculate total density
%     #   2. loop through cells, calculate total required number of macro particles
%     #       2.1. calculate the average density in that cell (i just took the center density as rough estimate)
%     #       2.2. calculate the cell volume (area) dx*dy
%     #       2.3. calculate how large part of the total density should go into that cell
%     #       2.4. calculate how many of the total amount of particles this corresponds to        
%     #   3. loop through cells, 
%     #       3.1. initiate the right number of particles in each cell
%     #       3.2. calculate weight: macro-particle weight = (species density x cell hypervolume)/number of macro-particles per cell  
% 
%     print('initializing particles')
%     n_center = np.empty([int(nx)-1, int(ny)-1])
%     for i in np.arange(0,int(nx)-1):
%         for j in np.arange(0,int(ny)-1):
%             n_center[i,j] = nfun(x_center[i],y_center[j])
% 
%     Ntot = sum(n_center.flatten())*dx*dy # approximate total number of micro particles
%     # create empty arrays
%     N         = np.empty([int(nx)-1, int(ny)-1]) # N - number of micro particles in each cell
%     Nfrac     = np.empty([int(nx)-1, int(ny)-1]) # N/Ntot - fraction of total number of microparticles in each cell
%     M         = np.empty([int(nx)-1, int(ny)-1]) # M - number of macro particles in each cell
% 
%     # First calculate the required number of macro particles M
%     for i in np.arange(0,int(nx)-1):
%         for j in np.arange(0,int(ny)-1):
%             # Calculate how many macro particles should go in each cell
%             cell_density = nfun(x_center[i],y_center[j])
%             cell_hypervolume = dx*dy
%             N[i,j] = cell_density*cell_hypervolume # density x area, number of micro particles in this cell, non-integer
%             Nfrac[i,j] = N[i,j]/Ntot # fraction of micro particles that is in this cell, i.e. for a uniform density, all numbers are the same
%             M[i,j] = Nfrac[i,j]*Mtot_target # number of macro particles in each cell, non-integer
%             
%     Mceil = np.ceil(M) # rounded up to have entire particles
%     Mceil = Mceil.astype(int) # make into int
%     Mtot = int(sum(Mceil.flatten())) # total number of macro particles
%     print('Mtot = ' +  str(Mtot))
% 
%     # Now initialize the arrays with the right dimensions Mtot and run loop again filling upp the arrays        
%     x_part    = np.zeros(Mtot)
%     y_part    = np.zeros(Mtot)
%     weight    = np.zeros(Mtot)
%
%     p_count = 0
%     for i in np.arange(0,int(nx)-1):
%         #print(i)
%         for j in np.arange(0,int(ny)-1):
%             p_count = p_count + Mceil[i,j] # increment particle counter
%             # Calculate weigths of particles, microparticles/macroparticles
%             # weight_cell = N[i,j]/Mceil[i,j]
% 
%             # Calculate particle positions      
%             xlo = x[i]
%             xhi = x[i+1]
%             ylo = y[j]
%             yhi = y[j+1]                                      
% 
%             istart = p_count - Mceil[i,j]
%             istop = p_count
%
%             x_part[istart:istop] = np.random.uniform(xlo,xhi,size=Mceil[i,j])
%             y_part[istart:istop] = np.random.uniform(ylo,yhi,size=Mceil[i,j])
%             weight[istart:istop] = N[i,j]/Mceil[i,j]
% 
%     Nfinal = p_count
%     x_y_weight = np.empty([3,Nfinal])
%     x_y_weight[0,:] = x_part
%     x_y_weight[1,:] = y_part
%     x_y_weight[2,:] = weight
% 
%     return x_y_weight
% ----------------------------------------------------------------------