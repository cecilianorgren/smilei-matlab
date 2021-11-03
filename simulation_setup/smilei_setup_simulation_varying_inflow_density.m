% Overall simulation stuff
mime = 25;
resx = 4;
resy = 4;
cell_length = [1/resx 1/resy];

box_size = [448 56]*sqrt(mime);
%box_size = [56 56]*sqrt(mime);
x0 = box_size(1)/2;
y0 = box_size(2)/2;
nxny = box_size./cell_length;
x_grid = 0:1/resx:box_size(1);
y_grid = 0:1/resy:box_size(2);


% Harris sheet 
n0 = 1;
nbg_harris = 0.001;
L = 2*sqrt(mime);
fn_Harris = @(x,y) nbg_harris + n0./cosh((y-y0)/L).^2; 

% Use symbolic expression to construct density variations (adding different
% functions)
syms xvar yvar
lin = 0.5*sqrt(mime); % transition length scale

% Cold particles from the top
n_cold_bg = 0.05;
n_cold_var = 0.1; % amplitude of transition
% locations of transitions
y_trans = (3:3:24)*sqrt(mime);
f_ctop_all = 0;
for itrans = 1:numel(y_trans)
  if mod(itrans,2) sign_f = 1; else sign_f = -1; end
  f_tmp = sign_f*0.5*n_cold_var*(1 + tanh(((yvar-y0)-y_trans(itrans))/lin));
  %f_tmp = 0.5*n_cold_var*(1 + sin((abs(yvar-y0))/2));
  f_ctop{itrans} = f_tmp;
  f_ctop_all = f_ctop_all + f_tmp;
end
fn_cold_var_top = matlabFunction(f_ctop_all);
fn_cold_bg = 0.5*n_cold_bg*(1+tanh((yvar-y0)/L));
fn_top = matlabFunction(f_ctop_all + fn_cold_bg,'vars',[xvar yvar]);

% Cold particles from the top
lin = 0.5*sqrt(mime); % transition length scale
n_cold_bg = 0.05;
n_cold_var = 0.1; % amplitude of transition
% locations of transitions
f_cbot_all = 0;
for itrans = 1:numel(y_trans)
  if mod(itrans,2) sign_f = 1; else sign_f = -1; end
  f_tmp = sign_f*0.5*n_cold_var*(1 + tanh((-(yvar-y0)-y_trans(itrans))/lin));
  %f_tmp = 0.5*n_cold_var*(1 + sin((abs(yvar-y0))/2));
  f_ctop{itrans} = f_tmp;
  f_cbot_all = f_cbot_all + f_tmp;
end
fn_cold_var_bot = matlabFunction(f_cbot_all);
fn_cold_bg = 0.5*n_cold_bg*(1+tanh(-(yvar-y0)/L));
fn_bot = matlabFunction(f_cbot_all + fn_cold_bg,'vars',[xvar yvar]);

% About how many particles we want
%M_target = [5e6,5e6,5e8,5e8,5e8,5e8];

M_target_harris = 2e7;
M_target_cold = 2e8;
M_target = [M_target_harris,M_target_harris,M_target_cold,M_target_cold,M_target_cold,M_target_cold];

M_tot_target = sum(M_target);

% Collect distributions for all species
n_dist = {fn_Harris,fn_Harris,fn_top,fn_top,fn_bot,fn_bot};

% Just do two species to speed up
%M_target_harris = 1e7;
%M_target_cold = 1e8;
%M_target = [M_target_harris,M_target_cold];
%M_tot_target = sum(M_target);
%n_dist = {fn_Harris,fn_top};

% Initalize particles
disp('Initializing particles.')
particles  = smilei_initialize_particles_write_h5(M_target,n_dist,[0,box_size(1),0,box_size(2)],nxny);
disp('Ready.')
%% Save to file
% On betzy, filepath should be something like:
% particle_position.h5/particles/species1
% particle_position.h5/particles/species2
% ...
%
% To move file to betzy scp
% scp /Users/cno062/Data/SMILEI/initialize_particles/particle_position.h5 cno062@betzy.sigma2.no:"/cluster/projects/nn9496k/Smilei/cn_cold_dipolarization/" 
disp('Writing h5 file.')
filename = '/Users/cno062/Data/SMILEI/initialize_particles/particle_position.h5';
%filename = '/Users/cno062/Data/SMILEI/initialize_particles/particle_position_test.h5';
delete(filename);
for iSpecies = 1:numel(particles)
  fprintf('\n Species %g',iSpecies)
  group = sprintf('/particles/species%g/',iSpecies);
  data = particles(iSpecies).weight;
  dataset = [group 'weight'];
  h5create(filename,dataset,numel(data)) % using size(data) is not compatible with Smilei, because then the data is "2D, and not 1D as required"
  h5write(filename,dataset,data)
  for comp = ['x','y','z']
    data = particles(iSpecies).(comp);
    dataset = [group 'position/' comp];
    h5create(filename,dataset,numel(data))
    h5write(filename,dataset,data)
  end
end
disp('Done writing h5 file.')

% Debugging from before
% % Check data_space
% filename = '/Users/cno062/Data/SMILEI/initialize_particles/particle_position.h5';
% fid = H5F.open(filename,'H5P_DEFAULT');
% dset_id = H5D.open(fid,'/particles/species1/weight');
% space = H5D.get_space(dset_id);
% [~,dims] = H5S.get_simple_extent_dims(space);
% H5S.close(space);
% H5D.close(dset_id);
% H5F.close(fid);
% %%
% filename = '/Users/cno062/Data/SMILEI/initialize_particles/particle_position.h5';
% fid = H5F.open(filename);
% dset_id = H5D.open(fid,'/particles/species1/weight');
% space_id = H5D.get_space(dset_id);
% [ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);
% matlab_dims = fliplr(h5_dims);
%% Plot results
nrows = 2;%numel(M_target);
ncols = 2;
h = setup_subplots(nrows,ncols);
isub = 1;

for iSpecies = 2:3%1:numel(M_target)
  if 0 % Individual particle locations  
    hca = h(isub); isub = isub + 1;
    plot(hca,particles(iSpecies).x,particles(iSpecies).y,'.')  
  end
  if 1 % Particles per bin
    hca = h(isub); isub = isub + 1;
    [count,edges,mid,loc] = histcn([particles(iSpecies).x,particles(iSpecies).y],x_grid,y_grid);
    pcolor(hca,mid{1:2},count')
    shading(hca,'flat');
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = {'macro particles per cell',sprintf('total = %.2f million',sum(count(:))*1e-6)};
  end
  if 1 % Individual particles binned
    hca = h(isub); isub = isub + 1;
    [count,edges,mid,loc] = histcn([particles(iSpecies).x,particles(iSpecies).y],x_grid,y_grid,'AccumData',particles(iSpecies).weight);
    pcolor(hca,mid{1:2},count'/prod(cell_length))
    shading(hca,'flat');
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'Binned particle weight = density' ;
  end
end

hlinks = linkprop(h,{'XLim','YLim'});
h(1).XLim = [0,box_size(1)];
h(1).YLim = [0,box_size(2)];