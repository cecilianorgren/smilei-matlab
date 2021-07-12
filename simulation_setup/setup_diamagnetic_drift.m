%diamagnetic_drit
L = 1;
Lp = 0.2;
zp = 4:3:14;
B0 = 1;
n0 = 1;
nb = 0.3;
T = 0.5;
Ptot = (2*B0^2)/2;
npert = @(z) nb*(1 + tanh((z-zp(1))/Lp) - tanh((z-zp(2))/Lp) + tanh((z-zp(3))/Lp) - tanh((z-zp(4))/Lp));
n = @(z) n0*cosh(z/L).^(-2) + npert(z);
Pp = @(z) n(z)*T;

Pppert = @(z) npert(z)*T;
PBpert = @(z) Ptot - Pppert(z);
%Bmagpert = @(z) sqrt(2*PBpert(z));
Bxpert = @(z) sqrt(PBpert(z));
Bypert = @(z) sqrt(PBpert(z));

%Bpert0 = sqrt((Ptot - T*nb)*2);
%Bpert = @(z) sqrt((Ptot - T*nb*npert(z))*2)/Bpert0;
Bx = @(z) B0*tanh(z/L).*Bxpert(z);
By = @(z) B0*z./z.*Bypert(z);


PB = @(z) (Bx(z).^2 + By(z).^2)/2;

z = linspace(-15,15,1000);

shear = atand(Bx(z)./By(z));

h = setup_subplots(6,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,z,Bx(z),z,By(z))
legend(hca,'B_x','B_y')
hca.YLim = 1*[-1.5 1.5];
grid(hca,'on')

hca = h(isub); isub = isub + 1;
plot(hca,z,Bmagpert(z),z,Bxpert(z),z,Bypert(z))
legend(hca,'|B|^{pert}','B_{x}^{pert}','B_{y}^{pert}')
hca.YLim = [-1.5 1.5];
grid(hca,'on')

if 1
hca = h(isub); isub = isub + 1;
plot(hca,z,shear)
legend(hca,'\theta')
grid(hca,'on')
end

hca = h(isub); isub = isub + 1;
plot(hca,z,n(z))
legend(hca,'n')
grid(hca,'on')
hca.YLim(1) = 0;

hca = h(isub); isub = isub + 1;
plot(hca,z,Pp(z),z,PB(z),z,Pp(z)+PB(z))
legend(hca,'P_p','P_B','P_p+P_B','location','west')
%hca.YLim = 1*[-1.5 1.5];
grid(hca,'on')


hca = h(isub); isub = isub + 1;
plot(hca,z,Pp(z)./PB(z))
legend(hca,'\beta')
hca.YLim(1) = 0;
grid(hca,'on')

%%
theta = 90;
L = 1;
Lp = 0.2;
zp = 4:3:14;
%B0 = 1;
n0 = 1;
nb = 0.1;
T = 0.5;
npert = @(z) 1 + tanh((z-zp(1))/Lp) - tanh((z-zp(2))/Lp) + tanh((z-zp(3))/Lp) - tanh((z-zp(4))/Lp);
n = @(z) n0*cosh(z/L).^(-2) + nb*npert(z);

Pp = @(z) n(z)*T;
Pppert = @(z) npert(z)*T;

PBpert = @(z) Ptot - Pppert(z);

Bmagpert = @(z) sqrt(2*PBpert(z));

Bxpert = @(z) Bmagpert(z);
Bypert = @(z) Bmagpert(z);


%Bpert0 = sqrt((Ptot - T*nb)*2);
%Bpert = @(z) sqrt((Ptot - T*nb*npert(z))*2)/Bpert0;
%Bx = @(z) B0*tanh(z/L).*Bpert(z);
%By = @(z) B0*z./z.*Bpert(z);


%PB = @(z) (Bx(z).^2 + By(z).^2)/2;

z = linspace(-15,15,1000);

%shear = atand(Bx(z)./By(z));

h = setup_subplots(6,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,z,n(z))
legend(hca,'n')
grid(hca,'on')
hca.YLim(1) = 0;

hca = h(isub); isub = isub + 1;
plot(hca,z,Pp(z),z,PB(z),z,Pp(z)+PB(z))
legend(hca,'P_p','P_B','P_p+P_B','location','west')
hca.YLim = [0 1.5];
grid(hca,'on')

hca = h(isub); isub = isub + 1;
plot(hca,z,Bmag(z),z,Bx(z),z,By(z))
legend(hca,'|B|','B_x','B_y')
hca.YLim = 1*[-1.5 1.5];
grid(hca,'on')

%%

hca = h(isub); isub = isub + 1;
plot(hca,z,Bx(z),z,By(z))
legend(hca,'B_x','B_y')
hca.YLim = 1*[-1.5 1.5];
grid(hca,'on')

hca = h(isub); isub = isub + 1;
plot(hca,z,Bpert(z))
legend(hca,'B_{pert}')
%hca.YLim = 1*[-1.5 1.5];
grid(hca,'on')

if 1
hca = h(isub); isub = isub + 1;
plot(hca,z,shear)
legend(hca,'\theta')
grid(hca,'on')
end


hca = h(isub); isub = isub + 1;
plot(hca,z,Pp(z),z,PB(z),z,Pp(z)+PB(z))
legend(hca,'P_p','P_B','P_p+P_B','location','west')
%hca.YLim = 1*[-1.5 1.5];
grid(hca,'on')


hca = h(isub); isub = isub + 1;
plot(hca,z,Pp(z)./PB(z))
legend(hca,'\beta')
hca.YLim(1) = 0;
grid(hca,'on')