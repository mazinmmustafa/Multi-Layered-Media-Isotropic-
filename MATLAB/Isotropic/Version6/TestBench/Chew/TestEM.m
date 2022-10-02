close all; clear; clc;
%%
Ns          =   1000;
[m,~,~,~,~,~,~,~,~]=Units();
%% Definitions
theta0      =   20;
phi0        =   30;
z_          =   -1.4*m;
y           =   1*m;
z           =   -0.3*m;
x_min    	=   -3*m;
x_max    	=   3*m;
%% Corrections
x           =   linspace(x_min,x_max,Ns);
y           =   y*ones(size(x));
[phi,rho]   =   cart2pol(x,y);
%% Dipole Components
[Mx,My,Mz]=M(theta0,phi0);
%% Compute Fields
Ex          =   zeros(1,Ns);
Ey          =   zeros(1,Ns);
Ez          =   zeros(1,Ns);
tic;
for i=1:Ns
    Ex(1,i)     =   GEMxx(rho(i),phi(i),z,z_)*Mx+GEMxy(rho(i),phi(i),z,z_)*My+GEMxz(rho(i),phi(i),z,z_)*Mz;
    Ey(1,i)     =   GEMyx(rho(i),phi(i),z,z_)*Mx+GEMyy(rho(i),phi(i),z,z_)*My+GEMyz(rho(i),phi(i),z,z_)*Mz;
    Ez(1,i)     =   GEMzx(rho(i),phi(i),z,z_)*Mx+GEMzy(rho(i),phi(i),z,z_)*My+0*GEMzz(rho(i),phi(i),z,z_)*Mz;
    clc;
    fprintf('Step:\t%0.0f/100\n',i*100/Ns);
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%%
figure()
plot(x/m,abs(Ex),'k','LineWidth',1)
xlabel('$x$ [m]','Interpret','Latex')
ylabel('$|E_{x}|$ [V/m]','Interpret','Latex')
title('Electric Field','Interpret','Latex')
xlim([x_min x_max]/m)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ex' -pdf -transparent
%%
figure()
plot(x/m,abs(Ey),'k','LineWidth',1)
xlabel('$x$ [m]','Interpret','Latex')
ylabel('$|E_{y}|$ [V/m]','Interpret','Latex')
title('Electric Field','Interpret','Latex')
xlim([x_min x_max]/m)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ey' -pdf -transparent
%%
figure()
plot(x/m,abs(Ez),'k','LineWidth',1)
xlabel('$x$ [m]','Interpret','Latex')
ylabel('$|E_{z}|$ [V/m]','Interpret','Latex')
title('Electric Field','Interpret','Latex')
xlim([x_min x_max]/m)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ez' -pdf -transparent
%% Dipole Components
function[Mx,My,Mz]=M(theta0,phi0)
theta0      =   deg2rad(theta0);
phi0        =   deg2rad(phi0);
Mx          =   sin(theta0)*cos(phi0);
My          =   sin(theta0)*sin(phi0);
Mz          =   cos(theta0);
end
%%