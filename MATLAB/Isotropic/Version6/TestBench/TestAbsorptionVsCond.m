close all; clear; clc;
%%
global var
%%
Ns          =   1e3;
mS          =   1e-3;
%%
sigma_max   =   100*mS;
var         =   0;
theta       =   0;
[N,~,~,~,~,~,Data]=Configs();
k1          =   Data(2,4);
kN          =   Data(N+1,4);
sigma       =   linspace(0,sigma_max,Ns);
%% Gamma Down
kpi         =   k1*sin(theta*pi/180);
tau_e       =   zeros(1,Ns);
tau_h       =   zeros(1,Ns);
Gamma_e     =   zeros(1,Ns);
Gamma_h     =   zeros(1,Ns);
for i=1:Ns
    var      	=   sigma(i);
    [Gamma_e(1,i),Gamma_h(1,i)]=Refl(kpi,'D',1,1);
    [tau_e(1,i),tau_h(1,i)]=tau_(kpi,1,1,'D',1,N,1);
end
[~,~,Ze1,Zh1]=kpParameters(kpi,1,1);
[~,~,ZeN,ZhN]=kpParameters(kpi,N,1);
T_e         =   real(Ze1./ZeN).*abs(tau_e).^2;
T_h         =   real(Zh1./ZhN).*abs(tau_h).^2;
R_e         =   abs(Gamma_e).^2;
R_h         =   abs(Gamma_h).^2;
A_e         =   1-T_e-R_e;
A_h         =   1-T_h-R_h;
figure()
hold on
plot(sigma/mS,A_e,'k','LineWidth',1)
plot(sigma/mS,A_h,'--k','LineWidth',1)
hold off
xlabel('$\sigma_s$ [mS]','Interpret','Latex')
ylabel('$|A|^2$','Interpret','Latex')
title('Absorption','Interpret','Latex')
legend('$e$','$h$','Interpreter','Latex','Location','NorthEast')
xlim([0 sigma_max/mS])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\A_sigma' -pdf -transparent
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

