close all; clear; clc;
%%
Ns          =   1e3;
%%
[~,~,~,um,~,~,~,~,~]=Units();
[N,~,k0,~,~,~,Data]=Configs();
kp          =   0.7*k0;
z_          =   500*um;
z0        	=   Data(1,1);
zN        	=   Data(N+1,1);
z           =   linspace(zN,z0,Ns);
%% 
V_e         =   zeros(1,Ns);
I_e         =   zeros(1,Ns);
V_h         =   zeros(1,Ns);
I_h         =   zeros(1,Ns);
for i=1:Ns
	[V_e(1,i),I_e(1,i),V_h(1,i),I_h(1,i)]=TLGF(kp,'i',z(i),z_,1);   
end
%% Voltage e
figure()
hold on
plot(z/um,real(V_e),'k','LineWidth',1)
plot(z/um,imag(V_e),'--k','LineWidth',1)
y_range     =   PlotLayers(um);
hold off
xlabel('$z$ [$\mu$m]','Interpret','Latex')
ylabel('$V^{e}_{v}$','Interpret','Latex')
title('Voltage','Interpret','Latex')
legend('Re','Im','Interpreter','Latex')
xlim([zN z0]/um)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ve' -pdf -transparent
%% Current e 
figure()
hold on
plot(z/um,real(I_e),'k','LineWidth',1)
plot(z/um,imag(I_e),'--k','LineWidth',1)
y_range     =   PlotLayers(um);
hold off
xlabel('$z$ [$\mu$m]','Interpret','Latex')
ylabel('$I^{e}_{v}$','Interpret','Latex')
title('Current','Interpret','Latex')
legend('Re','Im','Interpreter','Latex')
xlim([zN z0]/um)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ie' -pdf -transparent
%% Voltage h
figure()
hold on
plot(z/um,real(V_h),'k','LineWidth',1)
plot(z/um,imag(V_h),'--k','LineWidth',1)
y_range     =   PlotLayers(um);
hold off
xlabel('$z$ [$\mu$m]','Interpret','Latex')
ylabel('$V^{h}_{v}$','Interpret','Latex')
title('Voltage','Interpret','Latex')
legend('Re','Im','Interpreter','Latex')
xlim([zN z0]/um)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Vh' -pdf -transparent
%% Current h 
figure()
hold on
plot(z/um,real(I_h),'k','LineWidth',1)
plot(z/um,imag(I_h),'--k','LineWidth',1)
y_range     =   PlotLayers(um);
hold off
xlabel('$z$ [$\mu$m]','Interpret','Latex')
ylabel('$I^{h}_{v}$','Interpret','Latex')
title('Current','Interpret','Latex')
legend('Re','Im','Interpreter','Latex')
xlim([zN z0]/um)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ih' -pdf -transparent
%%