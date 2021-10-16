%% Development
dirpath = '/Users/cecilia/Discs/betzy/Smilei/cold_dipolarization_CN/cold_momentum/';
filename = 'ParticleBinning77.h5';

info = h5info([dirpath filename]);
name = info.Attributes(2).Value;


data = h5read([dirpath filename],'/timestep00001300');


%% Plot particle binning results
nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

iter = 1400;


if 0 % id = 36:41, ekin_dist_ 
  timestep = sprintf('/timestep%08.0f',iter);
  hca = h(isub); isub = isub + 1;
  id_all = 36:41;
  species = cell(6,1);
  %plot(hca,nan,nan)
  hold(hca,'on')
  for ii = 1:numel(id_all)
    id = id_all(ii);
    filename = sprintf('ParticleBinning%g.h5',id);
    info = h5info([dirpath filename]);
    data = h5read([dirpath filename],timestep);
    name = h5readatt([dirpath filename],'/','name');
    dep1_str = h5readatt([dirpath filename],'/','axis0');
    dep1_tok = split(dep1_str);
    if str2num(dep1_tok{5}) == 1
      dep1 = logspace(str2num(dep1_tok{2}),str2num(dep1_tok{3}),str2num(dep1_tok{4}));
    else
      dep1 = linspace(str2num(dep1_tok{2}),str2num(dep1_tok{3}),str2num(dep1_tok{4}));
    end
    species{ii} = h5readatt([dirpath filename],'/','species');
    hl(ii) = plot(hca,dep1,squeeze(data));
    hca.XLabel.String = dep1_tok{1};
    hca.YLabel.String = name(1:9);
    hca.YLabel.Interpreter = 'none';
    if ii == 1
      hold(hca,'on')
    elseif ii == numel(id_all)
      hold(hca,'off')
    end
  end
  hold(hca,'off')
  hca.XScale = 'log';
  hca.YScale = 'log';
  hleg = legend(hca,species,'location','best');
  hleg.Title.String = 'species';
end

if 1 % id = 54, vdf_px_ 
  hca = h(isub); isub = isub + 1; 
  id = 54+1*6+2+0;
  timestep = sprintf('/timestep%08.0f',iter);
  doLog = 1;
  filename = sprintf('ParticleBinning%g.h5',id);
  info = h5info([dirpath filename]);
  data = h5read([dirpath filename],timestep);
  name = h5readatt([dirpath filename],'/','name');
  dep0 = h5readatt([dirpath filename],'/','axis0');
  dep1 = h5readatt([dirpath filename],'/','axis1');
  dep2 = h5readatt([dirpath filename],'/','axis2');
  [dep0_str,dep0_val] = make_dep(dep0);
  [dep1_str,dep1_val] = make_dep(dep1);
  [dep2_str,dep2_val] = make_dep(dep2);
  
  if doLog
    imagesc(hca,dep0_val,dep2_val,log10(squeeze(data)))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = ['log10 ' name];
  else
    imagesc(hca,dep0_val,dep2_val,squeeze(data))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = name;
  end
  hca.XLabel.String = dep0_str;
  hca.YLabel.String = dep2_str;
  hcb.YLabel.Interpreter = 'none';
  hca.YDir = 'normal';
end
if 0 % id = 72-27, vdf_ekin_ 
  hca = h(isub); isub = isub + 1; 
  id = 72+2;
  timestep = sprintf('/timestep%08.0f',iter);
  doLog = 1;
  filename = sprintf('ParticleBinning%g.h5',id);
  info = h5info([dirpath filename]);
  data = h5read([dirpath filename],timestep);
  name = h5readatt([dirpath filename],'/','name');
  dep0 = h5readatt([dirpath filename],'/','axis0');
  dep1 = h5readatt([dirpath filename],'/','axis1');
  dep2 = h5readatt([dirpath filename],'/','axis2');
  [dep0_str,dep0_val] = make_dep(dep0);
  [dep1_str,dep1_val] = make_dep(dep1);
  [dep2_str,dep2_val] = make_dep(dep2);
  
  if doLog
    imagesc(hca,dep0_val,dep2_val,log10(squeeze(data)))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = ['log10 ' name];
  else
    imagesc(hca,dep0_val,dep2_val,squeeze(data))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = name;
  end
  hca.XLabel.String = dep0_str;
  hca.YLabel.String = dep2_str;
  hcb.YLabel.Interpreter = 'none';
  hca.YDir = 'normal';
end

if 0 % id = 77
  hca = h(isub); isub = isub + 1;
  id = 77;
  filename = sprint('ParticleBinning%g.h5',id);
  info = h5info([dirpath filename]);
  data = h5read([dirpath filename],'/timestep00001300');
  name = h5readatt([dirpath filename],'/','name');
  dep1 = h5readatt([dirpath filename],'/','axis0');
  dep2 = h5readatt([dirpath filename],'/','axis1');
  dep3 = h5readatt([dirpath filename],'/','axis2');
  imagesc(hca,dep1,dep3,squeeze(data))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = name;
end

%% Temporal evolution of partice distributions
dirpath = '/Users/cecilia/Discs/betzy/Smilei/cold_dipolarization_CN/cold_momentum/';
doLog = 1;
h = setup_subplots(3,4);
times = 1000:100:2000;
nt = numel(times);
ids = 54 + [0:3 6:9 12:15];
nid = numel(ids);
for it = 1:nt  
  isub = 1;
  timestep = sprintf('/timestep%08.0f',times(it));
  for id = 1:nid
    filename = sprintf('ParticleBinning%g.h5',ids(id));
    hca = h(isub); isub = isub + 1; 
    
    %info = h5info([dirpath filename]);
    data = h5read([dirpath filename],timestep);
    name = h5readatt([dirpath filename],'/','name');
    dep0 = h5readatt([dirpath filename],'/','axis0');
    dep1 = h5readatt([dirpath filename],'/','axis1');
    dep2 = h5readatt([dirpath filename],'/','axis2');
    [dep0_str,dep0_val] = make_dep(dep0);
    [dep1_str,dep1_val] = make_dep(dep1);
    [dep2_str,dep2_val] = make_dep(dep2);

    if doLog
      imagesc(hca,dep0_val,dep2_val,log10(squeeze(data)))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = ['log10 ' name];
    else
      imagesc(hca,dep0_val,dep2_val,squeeze(data))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = name;
    end
    hca.XLabel.String = dep0_str;
    hca.YLabel.String = dep2_str;
    hca.Title.String = timestep(2:end);
    hcb.YLabel.Interpreter = 'none';
    hca.YDir = 'normal';
    colormap(hca,pic_colors('candy4'))
    drawnow
  end  
  cn.print(['vx_vy_vz_' timestep(2:end)])
end