close all; clear; clc;
%%
Ns          =   1000;
[~,~,~,~,nm,~,~,~,~]=Units();
%% Definitions
rho         =   633*nm;
phi         =   0;
z_min    	=   -2000*nm;
z_max    	=   2000*nm;
z_          =   500*nm;
%% Corrections
phi         =   deg2rad(phi);
z           =   linspace(z_min,z_max,Ns);
%% Compute DGF
Gzz         =   zeros(1,Ns);
Gxz         =   zeros(1,Ns);
Gxx         =   zeros(1,Ns);
Gzx         =   zeros(1,Ns);
tic;
for i=1:Ns
    Gxx(1,i)    =   GEJxx(rho,phi,z(i),z_);
    Gzx(1,i)    =   GEJzx(rho,phi,z(i),z_);
    Gxz(1,i)    =   GEJxz(rho,phi,z(i),z_);
    Gzz(1,i)    =   GEJzz(rho,phi,z(i),z_);
    clc;
    fprintf('Step:\t%0.0f/100\n',i*100/Ns);
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Plot GEJxx, GEJzx
figure()
hold on
plot(z/nm,log10(abs(Gxx)),'k')
plot(z/nm,log10(abs(Gzx)),'--k')
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{xx}$','$G^{EJ}_{zx}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Result1' -pdf -transparent
%% Plot GEJxz, GEJzz
figure()
hold on
plot(z/nm,20*log10(abs(Gxz)),'k')
plot(z/nm,20*log10(abs(Gzz)),'--k')
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$20\log_{10}|G^{EJ}|$','Interpret','Latex')
title('DGF','Interpret','Latex')
legend('$G^{EJ}_{xz}$','$G^{EJ}_{zz}$','Interpreter','Latex','Location','NorthWest')
xlim([z_min z_max]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Result2' -pdf -transparent
%%