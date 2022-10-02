close all; clear; clc;
%%
Ns          =   1000;
[~,~,~,~,nm,~,~,~,~]=Units();
%% Definitions
rho         =   633*nm;
phi         =   45;
z_min    	=   -1000*nm;
z_max    	=   1000*nm;
z_          =   750*nm;
%% Corrections
phi         =   deg2rad(phi);
z           =   linspace(z_min,z_max,Ns);
%% Compute DGF
Gxz         =   zeros(1,Ns);
Gxx         =   zeros(1,Ns);
Gzx         =   zeros(1,Ns);
Gyx         =   zeros(1,Ns);
tic;
for i=1:Ns
    Gxx(1,i)    =   GEMxx(rho,phi,z(i),z_);
    Gzx(1,i)    =   GEMzx(rho,phi,z(i),z_);
    Gxz(1,i)    =   GEMxz(rho,phi,z(i),z_);
    Gyx(1,i)    =   GEMyx(rho,phi,z(i),z_);
    clc;
    fprintf('Step:\t%0.0f/100\n',i*100/Ns);
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Plot GEMxx, GEMzx
figure()
hold on
plot(z/nm,log10(abs(Gxx)),'k','LineWidth',1)
plot(z/nm,log10(abs(Gzx)),'--k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EM}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EM}_{xx}$','$G^{EM}_{zx}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Result1' -pdf -transparent
%% Plot GEMxz, GEMyx
figure()
hold on
plot(z/nm,log10(abs(Gxz)),'k','LineWidth',1)
plot(z/nm,log10(abs(Gyx)),'--k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EM}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EM}_{xz}$','$G^{EM}_{yx}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Result2' -pdf -transparent
%%