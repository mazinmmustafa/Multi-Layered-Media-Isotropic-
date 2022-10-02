close all; clear; clc;
%%
Ns          =   1000;
[~,~,~,~,nm,~,~,~,~]=Units();
%% Definitions
rho         =   633*nm;
phi         =   45;
z_min    	=   -1000*nm;
z_max    	=   +1000*nm;
z_          =   +750*nm;
%% Corrections
phi         =   deg2rad(phi);
z           =   linspace(z_min,z_max,Ns);
%% Compute DGF
Gxx         =   zeros(1,Ns);
Gxy         =   zeros(1,Ns);
Gxz         =   zeros(1,Ns);
Gyx         =   zeros(1,Ns);
Gyy         =   zeros(1,Ns);
Gyz         =   zeros(1,Ns);
Gzx         =   zeros(1,Ns);
Gzy         =   zeros(1,Ns);
Gzz         =   zeros(1,Ns);
tic;
for i=1:Ns
    Gxx(1,i)    =   GEJxx(rho,phi,z(i),z_);
    Gxy(1,i)    =   GEJxy(rho,phi,z(i),z_);
    Gxz(1,i)    =   GEJxz(rho,phi,z(i),z_);
    Gyx(1,i)    =   GEJyx(rho,phi,z(i),z_);
    Gyy(1,i)    =   GEJyy(rho,phi,z(i),z_);
    Gyz(1,i)    =   GEJyz(rho,phi,z(i),z_);
    Gzx(1,i)    =   GEJzx(rho,phi,z(i),z_);
    Gzy(1,i)    =   GEJzy(rho,phi,z(i),z_);
    Gzz(1,i)    =   GEJzz(rho,phi,z(i),z_);
    clc;
    fprintf('Step:\t%0.0f/100\n',i*100/Ns);
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Plot GEJxx
figure()
hold on
plot(z/nm,log10(abs(Gxx)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{xx}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Plot GEJxy
figure()
hold on
plot(z/nm,log10(abs(Gxy)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{xy}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Plot GEJxz
figure()
hold on
plot(z/nm,log10(abs(Gxz)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{xz}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Plot GEJyx
figure()
hold on
plot(z/nm,log10(abs(Gyx)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{yx}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Plot GEJyy
figure()
hold on
plot(z/nm,log10(abs(Gyy)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{yy}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Plot GEJyz
figure()
hold on
plot(z/nm,log10(abs(Gyz)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{yz}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Plot GEJzx
figure()
hold on
plot(z/nm,log10(abs(Gzx)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{zx}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Plot GEJzy
figure()
hold on
plot(z/nm,log10(abs(Gzy)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{zy}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Plot GEJzz
figure()
hold on
plot(z/nm,log10(abs(Gzz)),'k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{zz}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%%
