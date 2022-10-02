close all; clear; clc;
%%
global var
%%
Ns          =   2e2;
mS          =   1e-3;
%%
sigma_max   =   50*mS;
var         =   0;
[N,~,~,~,~,~,Data]=Configs();
k1          =   Data(2,4);
kN          =   Data(N+1,4);
theta       =   linspace(0,pi/2,Ns);
sigma       =   linspace(0,sigma_max,Ns);
[sigma,theta]=meshgrid(sigma,theta);
%% Gamma Down
kpi         =   k1*sin(theta);
tau_e       =   zeros(Ns,Ns);
tau_h       =   zeros(Ns,Ns);
Gamma_e     =   zeros(Ns,Ns);   
Gamma_h     =   zeros(Ns,Ns);
T_e         =   zeros(Ns,Ns);
T_h         =   zeros(Ns,Ns);
Ze1         =   zeros(Ns,Ns);
ZeN         =   zeros(Ns,Ns);
Zh1         =   zeros(Ns,Ns);
ZhN         =   zeros(Ns,Ns);
for i=1:Ns
    for ii=1:Ns
        var       	=   sigma(i,ii);
        [Gamma_e(i,ii),Gamma_h(i,ii)]=Refl(kpi(i,ii),'D',1,1);
        [tau_e(i,ii),tau_h(i,ii)]=tau_(kpi(i,ii),1,1,'D',1,N,1);
        [~,~,Ze1(i,ii),Zh1(i,ii)]=kpParameters(kpi(i,ii),1,1);
        [~,~,ZeN(i,ii),ZhN(i,ii)]=kpParameters(kpi(i,ii),N,1);
        T_e(i,ii)  	=   real(Ze1(i,ii)/ZeN(i,ii))*abs(tau_e(i,ii))^2;
        T_h(i,ii)  	=   real(Zh1(i,ii)/ZhN(i,ii))*abs(tau_h(i,ii))^2;
    end
end
R_e         =   abs(Gamma_e).^2;
R_h         =   abs(Gamma_h).^2;
A_e         =   1-T_e-R_e;
A_h         =   1-T_h-R_h;
%% Plot Reflectivity
figure()
% pcolor(sigma/mS,theta*180/pi,R_e)
xlim([0 sigma_max/mS])
ylim([0 90])
colormap jet
shading flat
colorbar
caxis([0 1])
xlabel('$\sigma_{s}$ [mS]','Interpret','Latex')
ylabel('$\theta_i$ [deg]','Interpret','Latex')
title('Reflectivity (TM)','Interpret','Latex')
pbaspect([1 1 1])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Data' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(sigma/mS,theta*180/pi,R_e,'Data_Matrix.dat');
%% Plot Transmittance
figure()
pcolor(sigma/mS,theta*180/pi,T_e)
xlim([0 sigma_max/mS])
ylim([0 90])
colormap jet
shading flat
colorbar
caxis([0 1])
xlabel('$\sigma_{s}$ [mS]','Interpret','Latex')
ylabel('$\theta_i$ [deg]','Interpret','Latex')
title('Transmittance (TM)','Interpret','Latex')
pbaspect([1 1 1])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Data' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(sigma/mS,theta*180/pi,T_e,'Data_Matrix.dat');
%% Plot Absorption
figure()
% pcolor(sigma/mS,theta*180/pi,A_e)
xlim([0 sigma_max/mS])
ylim([0 90])
colormap jet
shading flat
colorbar
caxis([0 1])
xlabel('$\sigma_{s}$ [mS]','Interpret','Latex')
ylabel('$\theta_i$ [deg]','Interpret','Latex')
title('Absorption (TM)','Interpret','Latex')
pbaspect([1 1 1])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Data' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(sigma/mS,theta*180/pi,A_e,'Data_Matrix.dat');
%% Plot Reflectivity
figure()
pcolor(sigma/mS,theta*180/pi,R_h)
xlim([0 sigma_max/mS])
ylim([0 90])
colormap jet
shading flat
colorbar
caxis([0 1])
xlabel('$\sigma_{s}$ [mS]','Interpret','Latex')
ylabel('$\theta_i$ [deg]','Interpret','Latex')
title('Reflectivity (TE)','Interpret','Latex')
pbaspect([1 1 1])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Data' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(sigma/mS,theta*180/pi,R_h,'Data_Matrix.dat');
%% Plot Transmittance
figure()
pcolor(sigma/mS,theta*180/pi,T_h)
xlim([0 sigma_max/mS])
ylim([0 90])
colormap jet
shading flat
colorbar
caxis([0 1])
xlabel('$\sigma_{s}$ [mS]','Interpret','Latex')
ylabel('$\theta_i$ [deg]','Interpret','Latex')
title('Transmittance (TE)','Interpret','Latex')
pbaspect([1 1 1])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Data' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(sigma/mS,theta*180/pi,T_h,'Data_Matrix.dat');
%% Plot Absorption
figure()
pcolor(sigma/mS,theta*180/pi,A_h)
xlim([0 sigma_max/mS])
ylim([0 90])
colormap jet
shading flat
colorbar
caxis([0 1])
xlabel('$\sigma_{s}$ [mS]','Interpret','Latex')
ylabel('$\theta_i$ [deg]','Interpret','Latex')
title('Absorption (TE)','Interpret','Latex')
pbaspect([1 1 1])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Data' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(sigma/mS,theta*180/pi,A_h,'Data_Matrix.dat');
%%
function[Vn_e,Vn_h]=tau_(kp,Vm_e,Vm_h,Dir,m,n,Sheet)
j           =   sqrt(-1);
%%
Vn_e       	=   Vm_e;
Vn_h       	=   Vm_h;
%% tau Up
if Dir=='U'
    for i=m-1:-1:n
        [~,Theta_n,~,~]=kpParameters(kp,i,Sheet);
        [~,Theta_m,~,~]=kpParameters(kp,i+1,Sheet);
        [Gamma_n_e,Gamma_n_h]=Refl(kp,'U',i,Sheet);
        [Gamma_m_e,Gamma_m_h]=Refl(kp,'U',i+1,Sheet);
        tau_e       =   (1+Gamma_m_e).*exp(-j*Theta_m)./(1+Gamma_n_e.*exp(-j*2*Theta_n));
        tau_h       =   (1+Gamma_m_h).*exp(-j*Theta_m)./(1+Gamma_n_h.*exp(-j*2*Theta_n));
        Vn_e    	=   tau_e.*Vn_e;
        Vn_h    	=   tau_h.*Vn_h;
    end
end
%% tau Down
if Dir=='D'
    for i=m+1:1:n
        [~,Theta_n,~,~]=kpParameters(kp,i,Sheet);
        [~,Theta_m,~,~]=kpParameters(kp,i-1,Sheet);
        [Gamma_n_e,Gamma_n_h]=Refl(kp,'D',i,Sheet);
        [Gamma_m_e,Gamma_m_h]=Refl(kp,'D',i-1,Sheet);
        tau_e       =   (1+Gamma_m_e).*exp(-j*Theta_m)./(1+Gamma_n_e.*exp(-j*2*Theta_n));
        tau_h       =   (1+Gamma_m_h).*exp(-j*Theta_m)./(1+Gamma_n_h.*exp(-j*2*Theta_n));
        Vn_e    	=   tau_e.*Vn_e;
        Vn_h    	=   tau_h.*Vn_h;
    end
end
end
%%

