close all; clear; clc;
%%
Ns          =   400;
[~,~,~,~,nm,~,~,~,~]=Units();
%% Definitions
theta0      =   0;
phi0        =   0;
z_          =   0*nm;
rho         =   633*nm;
phi         =   0;
z_min    	=   -2000*nm;
z_max    	=   2000*nm;
%% Corrections
phi         =   deg2rad(phi);
z           =   linspace(z_min,z_max,Ns);
%% Dipole Components
[Jx,Jy,Jz]=J(theta0,phi0);
%% Compute Fields
Ex          =   zeros(1,Ns);
Ey          =   zeros(1,Ns);
Ez          =   zeros(1,Ns);
tic;
for i=1:Ns
    Ex(1,i)     =   GEJxx(rho,phi,z(i),z_)*Jx+GEJxy(rho,phi,z(i),z_)*Jy+GEJxz(rho,phi,z(i),z_)*Jz;
    Ey(1,i)     =   GEJyx(rho,phi,z(i),z_)*Jx+GEJyy(rho,phi,z(i),z_)*Jy+GEJyz(rho,phi,z(i),z_)*Jz;
    Ez(1,i)     =   GEJzx(rho,phi,z(i),z_)*Jx+GEJzy(rho,phi,z(i),z_)*Jy+GEJzz(rho,phi,z(i),z_)*Jz;
    clc;
    fprintf('Step:\t%0.0f/100\n',i*100/Ns);
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Compute E Field
E           =   sqrt(abs(Ex).^2+abs(Ey).^2+abs(Ez).^2);
%% Plot E
figure()
hold on
plot(z/nm,20*log10(abs(E)),'k')
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$20\log_{10}|E|$','Interpret','Latex')
title('Gold Kretschmann','Interpret','Latex')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\E' -pdf -transparent
%% Dipole Components
function[Jx,Jy,Jz]=J(theta0,phi0)
theta0      =   deg2rad(theta0);
phi0        =   deg2rad(phi0);
Jx          =   sin(theta0)*cos(phi0);
Jy          =   sin(theta0)*sin(phi0);
Jz          =   cos(theta0);
end
%%