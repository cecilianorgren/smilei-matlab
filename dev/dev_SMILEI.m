% Implement 3D data reading
filepath = '/Users/cno062/tesla/software/Smilei4.5/Smilei/3d_turb/Fields0.h5';
namelist = '/Users/cno062/tesla/software/Smilei4.5/Smilei/3d_turb/Harris3d.py';
%sm = PIC(filepath,namelist);
sm = SMILEI(filepath,namelist);

%% Check if it corresponds to intended input
pic = sm(2);
tic
Bx = pic.Bx; toc
By = pic.By;
Bz = pic.Bz;
Ex = pic.Ex;
Ey = pic.Ey;
Ez = pic.Ez;
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