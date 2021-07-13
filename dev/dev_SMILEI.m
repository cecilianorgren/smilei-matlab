% Implement 3D data reading
filepath = '/Users/cno062/tesla/software/Smilei4.5/Smilei/3d_turb/Fields0.h5';
namelist = '/Users/cno062/tesla/software/Smilei4.5/Smilei/3d_turb/Harris3d.py';
%sm = PIC(filepath,namelist);
sm = SMILEI(filepath,namelist);

%% Check if it corresponds to intended input
pic = sm(sm.nt);
tic
Bx = pic.Bx; toc
By = pic.By;
Bz = pic.Bz;
%Ex = pic.Ex;
%Ey = pic.Ey;
%Ez = pic.Ez;
x = pic.xi;
y = pic.yi; 
z = pic.zi;
[X,Y,Z] = meshgrid(x,y,z);
Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);
%%
h = setup_subplots(1,1); % under pic-matlab
isub = 1;
% 3D density plot
hca = h(isub); isub = isub + 1;

isovalue = 0.5*max(Bmag(:));
surf1 = isosurface(X,Y,Z,permute(Bmag,[2 1 3]),isovalue);
p1 = patch(hca,surf1);
%isonormals(x,y,z,normalized_Free_Energy_map,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(hca,3); axis tight
camlight; lighting gouraud
isonormals(X,Y,Z,permute(Bmag,[2 1 3]),p1);
set(p1,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(hca,3); axis tight
camlight; lighting gouraud

hold(hca,'off')

%% Implement 2D particle_initialization
filepath = '/Users/cno062/betzy/stableSmilei/namelist_dev/Fields0.h5';
namelist = '/Users/cno062/betzy/stableSmilei/namelist_dev/Harris_dev.py';
diagPartBinning = '/Users/cno062/betzy/stableSmilei/namelist_dev/ParticleBinning*.h5';
%filepath = '/Users/cecilia/Discs/betzy/stableSmilei/namelist_dev/Fields0.h5';
%namelist = '/Users/cecilia/Discs/betzy/stableSmilei/namelist_dev/Harris_dev.py';
%sm = PIC(filepath,namelist);
sm = SMILEI(filepath,namelist,diagPartBinning);

%% Diagnostics
% We need something to run through the diagnostics files and investigate
% what they contain, because now the files just have numbers

% ekin_dist
info = h5info('/Users/cno062/betzy/stableSmilei/namelist_dev/ParticleBinning36.h5');
nd = numel(info.Datasets);
for id = 1:nd
  name = info.Datasets(id).Name;
  edf = h5read('/Users/cno062/betzy/stableSmilei/namelist_dev/ParticleBinning36.h5',['/' name]);
  plot(edf)
  title(name)
  pause(0.02)
end
%%
% vdf_x44_84_z0
h5path = '/Users/cno062/betzy/stableSmilei/namelist_dev/ParticleBinning37.h5';
info = h5info(h5path);
nd = numel(info.Datasets);
for id = 1:nd
  name = info.Datasets(id).Name;
  vdf = h5read(h5path,['/' name]);
  x = linspace(44,84,41);
  v = linspace(-1,1,100);
  fvx = squeeze(sum(sum(vdf,2),1));
  fvy = squeeze(sum(sum(vdf,3),1));
  fvz = squeeze(sum(sum(vdf,3),2));
  %
  h = setup_subplots(3,2,'vertical');
  isub = 1;
  
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,x,v,fvx); shading(hca,'flat'); hca.XLabel.String = 'v_x';
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,x,v,fvy); shading(hca,'flat'); hca.XLabel.String = 'v_y';
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,x,v,fvz); shading(hca,'flat'); hca.XLabel.String = 'v_z';
  
  hca = h(isub); isub = isub + 1;  
  plot(hca,v,squeeze(sum(fvx,2))); hca.XLabel.String = 'v_x';
  hca = h(isub); isub = isub + 1;  
  plot(hca,v,squeeze(sum(fvy,2))); hca.XLabel.String = 'v_y';
  hca = h(isub); isub = isub + 1;  
  plot(hca,v,squeeze(sum(fvz,2))); hca.XLabel.String = 'v_z';
  %plot(edf)
  %title(name)
  pause
end


 













