filePath = '/Users/cecilia/Data/PIC/data/Fields0.h5';
%fileInfo = h5info(filePath);

nPatches = fileInfo.Datasets.Dataspace;
iIter = h5readatt(filePath,'/','latest_timestep');


Bx = h5read(filePath,'/data/0000000000/Bx'); % de
pos_y = h5read(filePath,'/patch-000001/species-00-ion/Position-1'); % de


mom_x = h5read(filePath,'/patch-000000/species-00-ion/Momentum-0');
mom_y = h5read(filePath,'/patch-000000/species-00-ion/Momentum-1');
mom_z = h5read(filePath,'/patch-000000/species-00-ion/Momentum-2');

mom_xyz = [mom_x,mom_y,mom_z];

v_edges = linspace(-0.3,0.3,50);

[count,edges,mid,loc] = histcn(mom_xyz, v_edges, v_edges, v_edges);

imagesc(v_edges,v_edges,sum(count(:,:,:),3))


%%
pos_x = h5read(filePath,'/patch-000000/species-01-eon/Position-0'); % de
pos_y = h5read(filePath,'/patch-000000/species-01-eon/Position-1'); % de


mom_x = h5read(filePath,'/patch-000000/species-01-eon/Momentum-0');
mom_y = h5read(filePath,'/patch-000000/species-01-eon/Momentum-1');
mom_z = h5read(filePath,'/patch-000000/species-01-eon/Momentum-2');

mom_xyz = [mom_x,mom_y,mom_z];

v_edges = linspace(-0.3,0.3,50);

[count,edges,mid,loc] = histcn(mom_xyz, v_edges, v_edges, v_edges);

imagesc(v_edges,v_edges,sum(count(:,:,:),3))




