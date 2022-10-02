close all; clear; clc;
%%
Ns          =   1e4;
%%
[N,~,~,~,~,~,Data]=Configs();
k1          =   Data(2,4);
kN          =   Data(N+1,4);
theta       =   linspace(0,pi/2,Ns);
%% Gamma Down
kpi         =   k1*sin(theta);
[Gamma_e,Gamma_h]=Refl(kpi,'D',1,1);
figure()
hold on
plot(theta*180/pi,abs(Gamma_e).^2,'k','LineWidth',1)
plot(theta*180/pi,abs(Gamma_h).^2,'--k','LineWidth',1)
hold off
xlabel('$\theta_i$ [deg]','Interpret','Latex')
ylabel('$|\Gamma|^2$','Interpret','Latex')
title('Reflectivity','Interpret','Latex')
legend('$e$','$h$','Interpreter','Latex','Location','NorthWest')
ylim([0 1])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Gamma_D' -pdf -transparent
%% Gamma Up
kpi         =   kN*sin(theta);
[Gamma_e,Gamma_h]=Refl(kpi,'U',N,1);
figure()
hold on
plot(theta*180/pi,abs(Gamma_e).^2,'k','LineWidth',1)
plot(theta*180/pi,abs(Gamma_h).^2,'--k','LineWidth',1)
hold off
xlabel('$\theta_i$ [deg]','Interpret','Latex')
ylabel('$|\Gamma|^2$','Interpret','Latex')
title('Reflectivity','Interpret','Latex')
legend('$e$','$h$','Interpreter','Latex','Location','NorthWest')
ylim([0 1])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Gamma_U' -pdf -transparent
%%