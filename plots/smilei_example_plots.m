localuser = datastore('local','user');

filepath = '/Users/cno062/tesla/software/Smilei4.5/Smilei/cold_run/Fields0.h5';
namelist = '/Users/cno062/tesla/software/Smilei4.5/Smilei/cold_run/Harris_new.py';

filepath = '/Users/cno062/tesla/software/smilei/Smilei/baselinerun/newpert/Fields0.h5';
namelist = '/Users/cno062/tesla/software/smilei/Smilei/baselinerun/newpert/Harris_newpert.py';

filepath = '/Users/cecilia/Discs/betzy/Smilei/cold_dipolarization_CN/Fields0.h5';
namelist = '/Users/cecilia/Discs/betzy/Smilei/cold_dipolarization_CN/Harris_wcold.py';

filepath = '/Users/cecilia/Discs/betzy/Smilei/cold_dipolarization_CN/cold_momentum/Fields0.h5';
namelist = '/Users/cecilia/Discs/betzy/Smilei/cold_dipolarization_CN/cold_momentum/Harris_wcold_extended_diagnostics.py';
particlebinningpath = '/Users/cecilia/Discs/betzy/Smilei/cold_dipolarization_CN/cold_momentum/';

filepath = ['/Users/' localuser '/Data/PIC/Smilei/cn_open_boundary/Fields0.h5'];
namelist = ['/Users/' localuser '/Data/PIC/Smilei/cn_open_boundary/Harris_open_boundary.py'];
particlebinningpath = ['/Users/' localuser '/Data/PIC/Smilei/cn_open_boundary/'];

filepath = ['/Users/' localuser '/Data/SMILEI/cn_open_boundary/Fields0.h5'];
namelist = ['/Users/' localuser '/Data/SMILEI/cn_open_boundary/Harris_open_boundary.py'];
particlebinningpath = ['/Users/' localuser '/Data/SMILEI/cn_open_boundary/'];


filepath = '/Users/cno062/Discs/betzy2/stableSmilei/namelist_dev/Fields0.h5';
namelist = '/Users/cno062/Discs/betzy2/stableSmilei/namelist_dev/Harris_dev.py';

filepath = '/Users/cno062/Discs/betzy/Smilei/cn_cold_dipolarization/Fields0.h5';
namelist = '/Users/cno062/Discs/betzy/Smilei/cn_cold_dipolarization/Harris_varying_inflow_load_positions_from_h5.py';
particlebinningpath = '/Users/cno062/Discs/betzy/Smilei/cn_cold_dipolarization/';


filepath = '/Users/cno062/Discs/betzy/Smilei/cn_cold_dipolarization/data_tmp_apert0.75/Fields0.h5';
namelist = '/Users/cno062/Discs/betzy/Smilei/cn_cold_dipolarization/data_tmp_apert0.75/Harris_varying_inflow_load_positions_from_h5.py';
particlebinningpath = '/Users/cno062/Discs/betzy/Smilei/cn_cold_dipolarization/data_tmp_apert0.75/';

filepath = '/Users/cecilia/Discs/tesla/AGU/diamagnetic/Fields0.h5';
namelist = '/Users/cecilia/Discs/tesla/AGU/diamagnetic/Diamagnetic.py';
particlebinningpath = '/Users/cecilia/Discs/tesla/AGU/diamagnetic/';

filepath = '/Users/cecilia/Data/PIC/Smilei/AGU2/Fields0.h5';
namelist = '/Users/cecilia/Data/PIC/Smilei/AGU2/Diamagnetic.py';
particlebinningpath = '/Users/cecilia/Data/PIC/Smilei/AGU2/';

%sm = SMILEI(filepath,namelist,particlebinningpath);
sm = SMILEI(filepath,namelist,particlebinningpath);

%%
filepath = '/Users/cno062/tesla/software/Smilei4.5/Smilei/3d_turb/Fields0.h5';
namelist = '/Users/cno062/tesla/software/Smilei4.5/Smilei/3d_turb/Harris3d.py';
sm = PIC(filepath,namelist);


%%
pic = sm(1:2:sm.nt).twcilim([20 100]);
h = pic.movie({'ni'}','A',0.5,'cmap',{pic_colors('thermal')},'clim',{[0 2]});

%% Figure: Diagnostic 1
pic = sm(sm.nt);
xlim = pic.xi([1 end]);
zlim = pic.xi([1 end]);

varstrs = {'n(1)','n(3)','n(5)','Bx','Bz','Jy','vex','vix'}';
clims = {[0 1.2],[0 1.2],[0 1.2],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
cmap_th = pic_colors('thermal');
cmap_pst = pic_colors('pasteljet');
cmap_br = pic_colors('blue_red');

cmaps = {cmap_th,cmap_th,cmap_th,cmap_br,...
  cmap_br,cmap_br,cmap_br,cmap_br,cmap_br};

pic.xlim(xlim).zlim(zlim).plot_map(varstrs,'clim',clims,'cmap',cmaps)

%% Figure: Diagnostic 2
pic = sm(sm.nt);
xlim = pic.xi([1 end]);
zlim = pic.xi([1 end]);

varstrs = {'n(1)','n(3)','Ex','Ey','Ez'}';
clims = {[0 1.2],[0 1.2],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
cmap_th = pic_colors('thermal');
cmap_pst = pic_colors('pasteljet');
cmap_br = pic_colors('blue_red');

cmaps = {cmap_th,cmap_th,cmap_br,...
  cmap_br,cmap_br,cmap_br,cmap_br,cmap_br};

pic.xlim(xlim).zlim(zlim).plot_map(varstrs,'clim',clims,'cmap',cmaps)

%% Figure: Diagnostic 3, quick
pic = sm(sm.nt);
xlim = pic.xi([1 end]);
zlim = pic.xi([1 end]);

varstrs = {'n(1)','n(3)','Bx','Bz','Jy','vex','vix'}';
clims = {[0 1.2],[0 1.2],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
cmap_th = pic_colors('thermal');
cmap_pst = pic_colors('pasteljet');N
cmap_br = pic_colors('blue_red');

cmaps = {cmap_th,cmap_th,cmap_br,...
  cmap_br,cmap_br,cmap_br,cmap_br,cmap_br};

pic.xlim(xlim).zlim(zlim).plot_map(varstrs,'clim',clims,'cmap',cmaps)

%% Figure: Diagnostic 4
pic = sm(sm.nt);
xlim = pic.xi([1 end]);
zlim = pic.xi([1 end]);

varstrs = {'n(3)'}';
clims = {[0 1.2],[0 1.2],[0 1.2]};
cmap_th = pic_colors('thermal');
cmap_pst = pic_colors('pasteljet');
cmap_br = pic_colors('blue_red');

cmaps = {cmap_th,cmap_th,cmap_th};

pic.xlim(xlim).zlim(zlim).plot_map(varstrs,'clim',clims,'cmap',cmaps)

%% Figure: Cut of density
pic = sm(sm.nt);
xlim = pic.xi([20 21]);
zlim = pic.xi([1 end]);

varstrs = {{'n(1)','n(3)','n(5)'}}';
clims = {[0 1.2],[0 1.2],[0 1.2]};
cmap_th = pic_colors('thermal');
cmap_pst = pic_colors('pasteljet');
cmap_br = pic_colors('blue_red');

cmaps = {cmap_th,cmap_th,cmap_th};

pic.xlim(xlim).zlim(zlim).plot_line('y',varstrs)
