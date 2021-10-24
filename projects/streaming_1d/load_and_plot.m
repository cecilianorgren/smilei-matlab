%% Load object
directory = '/Users/cecilia/Smilei/streaming_instabilities/beam3/';
namelist = [directory 'tst1d_02_two_str_instability.py'];
directory = '/Users/cecilia/Smilei/streaming_instabilities/beam_buneman/';
namelist = [directory 'buneman.py'];

filepath = [directory 'Fields0.h5'];
particlebinningpath = directory;

%sm = SMILEI(filepath,namelist,particlebinningpath);

info = h5info(filepath);
nGroups = numel(info.Groups.Groups);
datasets = {info.Groups.Groups(1).Datasets.Name};
nDatasets = numel(datasets);
for iDataset = 1:nDatasets
  clear(datasets{iDataset});
end

for iGroup = 1:nGroups
  group = info.Groups.Groups(iGroup).Name;
  for iDataset = 1:nDatasets
    dataset = h5read(filepath,[group filesep datasets{iDataset}]);
    eval([datasets{iDataset} '(iGroup,:) = torow(dataset);' ])
  end
end

%% Diagonal particle binning
%% Temporal evolution of partice distributions
particlebinningpath = directory;
filename = sprintf('ParticleBinning%g.h5',1);
info_diag = h5info([particlebinningpath filesep filename]);
datasets = {info_diag.Datasets.Name};
nDatasets = numel(datasets);

clim = [0 3e-4];
doLog = 0;

h = setup_subplots(1,2);

for iDataset = 1:nDatasets
  dataset = info_diag.Datasets(iDataset).Name; 
  data = h5read([particlebinningpath filename],[filesep dataset]);
  if iDataset == 1
    data0 = data;
  end
  name = h5readatt([particlebinningpath filename],'/','name');
  dep0 = h5readatt([particlebinningpath filename],'/','axis0');
  dep1 = h5readatt([particlebinningpath filename],'/','axis1');  
  [dep0_str,dep0_val] = make_dep(dep0);
  [dep1_str,dep1_val] = make_dep(dep1);  

  isub = 1;
  if 1 % f(x,vx)
    hca = h(isub); isub = isub + 1;
    if doLog
      imagesc(hca,dep0_val,dep1_val,smooth2(log10(squeeze(data)),0))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = ['log10 ' name];
    else
      imagesc(hca,dep0_val,dep1_val,smooth2(squeeze(data),0))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = name;
    end
    %hca.CLim = clim;
    hca.XLabel.String = dep0_str;
  %    hca.YLabel.String = dep2_str;
    hca.Title.String = dataset;
    hcb.YLabel.Interpreter = 'none';
    hca.YDir = 'normal';
    colormap(hca,pic_colors('candy4'))
  end
  if 1 % <f(x)>
    hca = h(isub); isub = isub + 1;
    if doLog
      plot(dep1_val,mean(log10(squeeze(data)),2),dep1_val,mean(log10(squeeze(data0)),2))
      hca.YLabel.String = ['log10 ' name];    
    else
      plot(dep1_val,mean(squeeze(data),2),dep1_val,mean(squeeze(data0),2))
      hca.YLabel.String = name;
    end
    
    hca.XLabel.String = dep0_str;
    hca.Title.String = dataset;
    hca.YLabel.Interpreter = 'none';   
  end
  drawnow   
  pause(0.1)
  %cn.print(['vx_vy_vz_' timestep(2:end)])
end