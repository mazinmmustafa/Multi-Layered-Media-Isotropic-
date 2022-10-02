close all; clear; clc;
%%
global var
%%
Ns          =   6e2;
mS          =   1e-3;
%%
[~,lambda0,~,~,~,~,~]=Configs();
%% Definitions
x_min       =   -2*lambda0;   
x_max       =   2*lambda0;
z_min       =   -2*lambda0;
z_max       =   2*lambda0;
theta_i     =   45;
phi_i       =   180;
E_TM        =   1;
E_TE        =   0;
% var         =   5.305305305305*mS;
var         =   7.707707707708*mS;
%% Corrections
x           =   linspace(x_min,x_max,Ns);
y           =   0;
z           =   linspace(z_min,z_max,Ns);
[phi,rho]   =   cart2pol(x,y);
theta_i     =   deg2rad(theta_i);
phi_i       =   deg2rad(phi_i);
%% Compute Plane Wave
E_rho      	=   zeros(Ns,Ns);
E_phi      	=   zeros(Ns,Ns);
E_z      	=   zeros(Ns,Ns);
for i=1:Ns
    [E_rho(i,:),E_phi(i,:),E_z(i,:)]=EzCut(rho,phi,z(i),theta_i,phi_i,E_TM,E_TE);
end
%%
[x,z]       =   meshgrid(x,z);
E           =   sqrt(real(E_rho).^2+real(E_phi).^2+real(E_z).^2);
%% Plot Results
figure()
pcolor(x/lambda0,z/lambda0,E)
shading flat
colormap jet
colorbar
hold on
PlotLayersH(x_min,x_max,lambda0);
hold off
xlabel('$x/\lambda_0$','Interpret','Latex')
ylabel('$z/\lambda_0$','Interpret','Latex')
title('$|\textrm{Re}[E]|$','Interpret','Latex')
axis equal
axis([x_min x_max z_min z_max]/lambda0)
% caxis([0 7])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\PlaneWave' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(x/um,z/um,E,'Data_Matrix.dat');
%%

