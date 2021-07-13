function smilei_check_simulation(pic, varargin);


%% Check simulation for wierd stuff
if nargin == 1
    twcitime = pic.twci(end)
elseif nargin == 2
    [~,tind] = min(abs(pic.twci-varargin{1}'));
    twcitime = pic.twci(tind);
end
%%
figure; 

p = panel();
p.pack(3,2);
p(1,1).select(); hold on
imagesc(pic.xi,pic.yi,pic.twcilim(twcitime).ni'); set(gca,'YDir','normal'); hold on; set(gca,'layer','top'); colorbar
title(['dni   t= ' num2str(twcitime)])
%ylim([-25.2 25.2]);
%xlim([0 102]);
p(2,1).select(); hold on
imagesc(pic.xi,pic.yi,pic.twcilim(twcitime).Jz'); set(gca,'YDir','normal'); hold on; set(gca,'layer','top'); colorbar
title('jiz-jez')
%ylim([-25.2 25.2]);
%xlim([0 102]);
p(3,1).select(); hold on
imagesc(pic.xi,pic.yi,pic.twcilim(twcitime).Ez'); set(gca,'YDir','normal'); hold on; set(gca,'layer','top'); colorbar
title('ez')
%ylim([-25.2 25.2]);
%xlim([0 102]);
p(1,2).select(); hold on
imagesc(pic.xi,pic.yi,pic.twcilim(twcitime).Bz'); set(gca,'YDir','normal'); hold on; set(gca,'layer','top'); colorbar
title('by')
%ylim([-25.2 25.2]);
%xlim([0 102]);
p(2,2).select(); hold on
imagesc(pic.xi,pic.yi,pic.twcilim(twcitime).Bx'); set(gca,'YDir','normal'); hold on; set(gca,'layer','top'); colorbar
title('bx')
%ylim([-25.2 25.2]);
%xlim([0 102]);
p(3,2).select(); hold on
imagesc(pic.xi,pic.yi,pic.twcilim(twcitime).Ey'); set(gca,'YDir','normal'); hold on; set(gca,'layer','top'); colorbar
title('ey')
%ylim([-25.2 25.2]);
%xlim([0 102]);
%figuresize(35,24,'centimeters');
p.marginright = 22;
p.de.margin = 28;
% 
% figure(hh); hold on; plot(ze,dni(1600,:)); title(['dni   t= ' num2str(jj/wpewce/mass(1))])

%end