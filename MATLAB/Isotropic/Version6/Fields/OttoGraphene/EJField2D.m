close all; clear; clc;
%%
Ns          =   400;
[~,~,~,um,~,~,~,~,~]=Units();
%% Definitions
theta0      =   0;
phi0        =   0;
J           =   1;
x_          =   0*um;
y_          =   0*um;
z_          =   5*um;
y           =   0*um;
z_min     	=   -800*um;
z_max     	=   800*um;
x_min    	=   1000*um;
x_max    	=   2600*um;
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
E_abs    	=   sqrt(abs(Ex).^2+abs(Ey).^2+abs(Ez).^2);
E_real    	=   sqrt(real(Ex).^2+real(Ey).^2+real(Ez).^2);
%% Save Data
% DATA = real(Ex); save(strcat(pwd,'\Data\','ReEx.dat'),'DATA','-ascii');
% DATA = imag(Ex); save(strcat(pwd,'\Data\','ImEx.dat'),'DATA','-ascii');
% DATA = real(Ey); save(strcat(pwd,'\Data\','ReEy.dat'),'DATA','-ascii');
% DATA = imag(Ey); save(strcat(pwd,'\Data\','ImEy.dat'),'DATA','-ascii');
% DATA = real(Ez); save(strcat(pwd,'\Data\','ReEz.dat'),'DATA','-ascii');
% DATA = imag(Ez); save(strcat(pwd,'\Data\','ImEz.dat'),'DATA','-ascii');
%% Plot E
figure()
pcolor(x/um,z/um,20*log10(E_abs))
% hold on
% contour(x/um,z/um,20*log10(E_abs),30,'k')
% hold off
shading flat
colormap jet
colorbar
hold on
PlotLayersH(x_min,x_max,um);
hold off
xlabel('$x$ [$\mu$m]','Interpret','Latex')
ylabel('$z$ [$\mu$m]','Interpret','Latex')
title('$20\log_{10}|E|$','Interpret','Latex')
axis equal
axis([x_min x_max z_min z_max]/um)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\E' -pdf -transparent
%% Plot Ez
figure()
pcolor(x/um,z/um,20*log10(E_real))
% pcolor(x/um,z/um,20*log10(abs(real(Ez))))
shading flat
colormap jet
colorbar
hold on
PlotLayersH(x_min,x_max,um);
hold off
xlabel('$x$ [$\mu$m]','Interpret','Latex')
ylabel('$z$ [$\mu$m]','Interpret','Latex')
title('$20\log_{10}|\textrm{Re}[E]|$','Interpret','Latex')
axis equal
axis([x_min x_max z_min z_max]/um)
caxis([100 180])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\E' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(x/um,z/um,20*log10(E_real),'Data_Matrix.dat');
%%
function[Jx,Jy,Jz]=DipoleJ(theta0,phi0,J)
theta0      =   deg2rad(theta0);
phi0        =   deg2rad(phi0);
Jx          =   J*sin(theta0)*cos(phi0);
Jy          =   J*sin(theta0)*sin(phi0);
Jz          =   J*cos(theta0);
end
%%
