% Requriements:
% Number of patches must be a power of 2
% Cells/patch must be an even number
% Number of cells/number of patches = 3-10 is good

xi_goal = 150;         % about how big box you want in x/di
yi_goal = 50;          % about how big box you want in y/di
resx = 4;              % number of cells per de
resy = 4;              % number of cells per de
mime = 25;             % mass ratio
cells_per_patch = [6 8];  % how many cells per patch you want, changing this 
                       % will affect the final box size (xi_final, yi_final)
                       % What about different number of cells per patch in x
                       % and y directions?

xe_goal = xi_goal*sqrt(mime); % x/de
ye_goal = yi_goal*sqrt(mime); % y/de

n_cells_goal(1) = xe_goal*resx; % (#de)*(cells/de)
n_cells_goal(2) = ye_goal*resy; % (#de)*(cells/de)

% Find out what power of 2 gets you closest to cells_per_patch
power_all = 1:10;
cells_per_patch_goal = n_cells_goal'*2.^(-power_all);
[~,power_final] = min(abs(cells_per_patch_goal-repmat(reshape(cells_per_patch,[2 1]),1,size(cells_per_patch_goal,2))),[],2);
%[~,power_final] = find((cells_per_patch)>abs(cells_per_patch_goal),1,'first');
n_patches = 2.^power_all(power_final); % get number of patches

n_cells = n_patches.*cells_per_patch;

xe_final = n_cells(1)/resx;
ye_final = n_cells(2)/resy;

xi_final = xe_final/sqrt(mime);
yi_final = ye_final/sqrt(mime);

%disp(sprintf('n_cells = [%g, %g], cells_per_patch = [%g,%g], xyi_goal = [%g, %g], xyi_final = [%g, %g], n_patches = [%g, %g]',n_cells(1),n_cells(2),cells_per_patch(1),cells_per_patch(2),xi_goal,yi_goal,xi_final,yi_final,n_patches(1),n_patches(2)))
fprintf('---- SMILEI box and patches setup ---- \n%15s = [%6g,%6g] \n%15s = [%6g,%6g] \n%15s = [%6g,%6g] \n%15s = [%6g,%6g] \n%15s = [%6g,%6g] \n%15s = [%6g,%6g]\n','xyi_goal',xi_goal,yi_goal,'xye_goal',xe_goal,ye_goal,'cells_per_patch',cells_per_patch(1),cells_per_patch(2),'n_cells',n_cells(1),n_cells(2),'xyi_final',xi_final,yi_final,'n_patches',n_patches(1),n_patches(2))


