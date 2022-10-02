close all; clear; clc;
%%
Ns          =   400;
[~,~,~,~,nm,~,~,~,~]=Units();
%% Definitions
theta0      =   0;
phi0        =   0;
J           =   1;
x_          =   0*nm;
y_          =   0*nm;
z_          =   750*nm;
y           =   0*nm;
z_min     	=   -2000*nm;
z_max     	=   2000*nm;
x_min    	=   -2000*nm;
x_max    	=   2000*nm;
%% Corrections
[x,z]       =   meshgrid(linspace(x_min,x_max,Ns),linspace(z_min,z_max,Ns));
[phi,rho]   =   cart2pol(x-x_,y-y_);
%% Dipole Components
[Jx,Jy,Jz]=DipoleJ(theta0,phi0,J);
%% Compute Fields
Ex          =   zeros(Ns,Ns);
Ey          =   zeros(Ns,Ns);
Ez          =   zeros(Ns,Ns);
tic;
count       =   0;
for m=1:Ns
    for n=1:Ns
        Ex(m,n)     =   GEJxx(rho(m,n),phi(m,n),z(m,n),z_)*Jx+GEJxy(rho(m,n),phi(m,n),z(m,n),z_)*Jy+GEJxz(rho(m,n),phi(m,n),z(m,n),z_)*Jz;
        Ey(m,n)     =   GEJyx(rho(m,n),phi(m,n),z(m,n),z_)*Jx+GEJyy(rho(m,n),phi(m,n),z(m,n),z_)*Jy+GEJyz(rho(m,n),phi(m,n),z(m,n),z_)*Jz;
        Ez(m,n)     =   GEJzx(rho(m,n),phi(m,n),z(m,n),z_)*Jx+GEJzy(rho(m,n),phi(m,n),z(m,n),z_)*Jy+GEJzz(rho(m,n),phi(m,n),z(m,n),z_)*Jz;
        count       =   count+1;
        clc;
        fprintf('Step:\t%0.0f/100\n',count*100/Ns^2);
    end
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% E Field
E           =   sqrt(abs(Ex).^2+abs(Ey).^2+abs(Ez).^2);
%% Plot E
figure()
pcolor(x/nm,z/nm,20*log10(E))
shading flat
colormap jet
colorbar
hold on
PlotLayersH(x_min,x_max,nm);
hold off
xlabel('$x$ [nm]','Interpret','Latex')
ylabel('$y$ [nm]','Interpret','Latex')
title('$20\log_{10}|E|$','Interpret','Latex')
axis([x_min x_max z_min z_max]/nm)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\E' -pdf -transparent
%% Plot Ez
figure()
% pcolor(x/nm,z/nm,real(Ez))
pcolor(x/nm,z/nm,log10(abs(real(Ez))))
shading flat
colormap jet
colorbar
hold on
PlotLayersH(x_min,x_max,nm);
hold off
xlabel('$x$ [nm]','Interpret','Latex')
ylabel('$y$ [nm]','Interpret','Latex')
title('$\textrm{Re}[E_{z}]$','Interpret','Latex')
title('$\log_{10}|\textrm{Re}[E_{z}]|$','Interpret','Latex')
% axis([x_min x_max z_min z_max]/nm)
% caxis([-1e15 1e15])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ez' -pdf -transparent
%%
function[Jx,Jy,Jz]=DipoleJ(theta0,phi0,J)
theta0      =   deg2rad(theta0);
phi0        =   deg2rad(phi0);
Jx          =   J*sin(theta0)*cos(phi0);
Jy          =   J*sin(theta0)*sin(phi0);
Jz          =   J*cos(theta0);
end
%%
