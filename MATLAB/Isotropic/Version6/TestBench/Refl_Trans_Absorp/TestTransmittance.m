close all; clear; clc;
%%
Ns          =   1e3;
%%
[N,~,~,~,~,~,Data]=Configs();
k1          =   Data(2,4);
kN          =   Data(N+1,4);
theta       =   linspace(0,pi/2,Ns);
%% Gamma Down
kpi         =   k1*sin(theta);
tau_e       =   zeros(1,Ns);
tau_h       =   zeros(1,Ns);
[Gamma_e,Gamma_h]=Refl(kpi,'D',1,1);
for i=1:Ns
    [tau_e(1,i),tau_h(1,i)]=tau_(kpi(1,i),1,1,'D',1,N,1);
end
[~,~,Ze1,Zh1]=kpParameters(kpi,1,1);
[~,~,ZeN,ZhN]=kpParameters(kpi,N,1);
T_e         =   real(Ze1./ZeN).*abs(tau_e).^2;
T_h         =   real(Zh1./ZhN).*abs(tau_h).^2;
figure()
hold on
plot(theta*180/pi,T_e,'k','LineWidth',1)
plot(theta*180/pi,T_h,'--k','LineWidth',1)
hold off
xlabel('$\theta_i$ [deg]','Interpret','Latex')
ylabel('$|T|^2$','Interpret','Latex')
title('Transmittance','Interpret','Latex')
legend('$e$','$h$','Interpreter','Latex','Location','NorthWest')
xlim([0 90])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\tau_D' -pdf -transparent
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

