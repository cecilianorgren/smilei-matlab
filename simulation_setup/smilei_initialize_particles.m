function out  = smilei_initialize_particles(M_target_arr,n_dist,box_size,nxny)
% function out  = smilei_initialize_particles_write_h5(M_target,n_dist,box_size,[nx ny])
%   Calculates particle positions based on the total number of desired 
%   particles.
% 
% Returns table array with fields x, y, z, weight (one table for each species)
%   M_target - array with total number of target particles for each
%              species
%   n_dist   - cell array with distribution function
%   box_size - array with [xmin xmax ymin ymax]
%   [nx, ny] - number of cells in which to divide grid, does not have to be
%              the same as the final distribution, but preferably not
%              coarser

% Collect input
xmin = box_size(1); 
xmax = box_size(2); 
ymin = box_size(3); 
ymax = box_size(4);

nx = nxny(1);
ny = nxny(2);
  
nSpecies = numel(M_target_arr);

% Setup grid, for discretization of n
x = linspace(xmin,xmax,nx+1);
y = linspace(ymin,ymax,ny+1);
x0 = mean(x);
y0 = mean(y);
dx = x(2) - x(1);
dy = y(2) - y(1);
%cell_volume = dx*dy;
x_center = x(2:end)-0.5*dx;
y_center = y(2:end)-0.5*dy;
[X,Y] = ndgrid(x_center,y_center);

 
% Loop through species
fprintf('\n')
particles = [];
for iSpecies = 1:nSpecies
  tt = tic;
  fprintf(['Species %: ',iSpecies)
  % Density ditribution
  nfun = n_dist{iSpecies};
  
  % Target number of macro_particles int he simulation (for given species)
  M_tot_target = M_target_arr(iSpecies);

  % Discretize density distribution onto centers of grid.  
  N = nfun(X,Y);

  % Total number of microparticles in simulation (for given species)
  N_tot = sum(N(:)); 

  % Fraction of particles in each cell
  N_frac = N/N_tot;

  % Target number of macro particles in each cell
  M_target = M_tot_target*N_frac;

  % Final number of macro particles in each cell, round up M_target
  M = ceil(M_target);

  % Final particle weight, cellhypervolume*microparticles/macroparticles
  % https://smileipic.github.io/Smilei/units.html#weights
  W = dx*dy*N./M;

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
  
  particles(iSpecies).number_macro_particles = M_final;
  particles(iSpecies).x = part_x;
  particles(iSpecies).y = part_y;
  particles(iSpecies).z = part_y*0; % for 2D
  particles(iSpecies).weight = part_w;
  tt = toc(tt);
  fprintf(' Elapsed time is %g seconds\n',tt)
end

out = particles;

