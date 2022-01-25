%% Specify path
dirpath = '/Users/cecilia/Discs/betzy/Smilei/cn_cold_dipolarization/cold_momentum/';
dirpath = '/Users/cno062/Discs/betzy_base/work/users/paulten/diamagnetic/';
filename = 'TrackParticlesDisordered_eon.h5';

info = h5info([dirpath filename]);
name = info.Attributes(2).Value;

%data = h5read([dirpath filename],'/timestep00001300');

%% Read particles
%info.Groups(1).Groups(1).Groups(1).Groups(1).Groups(2).Datasets(1)
xlim = [0 1024];
ylim = [0 256];
nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

iter = 590;
species = 'eon_bgbumps';
filename = sprintf('TrackParticlesDisordered_%s.h5',species);
dataspace = info.Groups(1).Groups(1).Groups(1).Groups(1).Groups(2).Datasets(1).Dataspace.Size;
nLoadStart = dataspace;
nLoadStep = 1; nLoad = floor((dataspace-nLoadStart+1)/nLoadStep);

nLoadStart = 1;dataspace-60000;
nLoadStep = 1000; nLoad = floor((dataspace-nLoadStart+1)/nLoadStep);


if 1 % Postion (x,y)
  hca = h(isub); isub = isub + 1;
  iteration = sprintf('/timestep%08.0f',iter);
  dataset_x =  sprintf('/data/%010.0f/particles/%s/position/%s',iter,species,'x');
  dataset_y =  sprintf('/data/%010.0f/particles/%s/position/%s',iter,species,'y');
  data_x = h5read([dirpath filename],dataset_x,nLoadStart,nLoad,nLoadStep);
  data_y = h5read([dirpath filename],dataset_y,nLoadStart,nLoad,nLoadStep);
  plot(hca,data_x,data_y,'.')
  hca.XLabel.String = 'x/d_e';
  hca.YLabel.String = 'y/d_e';
  hca.Title.String = sprintf('Iteration = %010.0f, Species = %s, # particles plotted = %g',iter,species,numel(data_x));
  hca.XLim = xlim;
  hca.YLim = ylim;
end
