function varargout = hilbert_patches(nPatchesX,nPatchesY,xLim,yLim)
%%

%nOrderY = 

A = zeros(0,2);
B = zeros(0,2);
C = zeros(0,2);
D = zeros(0,2);

xlim = [-1 1];
ylim = [-1 1];
north = [ sum(xlim)/2  ylim(2)];
east  = [ xlim(2)  sum(ylim)/2];
south = [ sum(xlim)/2 ylim(1)];
west  = [ xlim(1)  sum(ylim)/2];

%north = [ 0  1];
%east  = [ 1  0];
%south = [ 0 -1];
%west  = [-1  0];

order = 6;
for n = 1:order
  AA = [B ; north ; A ; east  ; A ; south ; C];
  BB = [A ; east  ; B ; north ; B ; west  ; D];
  CC = [D ; west  ; C ; south ; C ; east  ; A];
  DD = [C ; south ; D ; west  ; D ; north ; B];

  A = AA;
  B = BB;
  C = CC;
  D = DD;
end

A = [0 0; cumsum(A)];
nPatches = numel(A);

hca = subplot(1,1,1);
plot(hca,A(:,1), A(:,2), 'clipping', 'off')
axis equal, axis off
%hold(hca,'on')
%%
for iPatch = 1:nPatches
  
  %plot(hca,[A(iPatch,1)+xlim(1), A(iPatch,2)+B(:,1), 'ro')
  
end
hold(hca,'off')
%axis equal, axis off