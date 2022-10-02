close all; clear; clc;
%%
Ns          =   1e3;
%%
[~,~,~,~,nm,~,~,~,~]=Units();
[N,~,k0,~,~,~,Data]=Configs();
kp          =   0.6*k0;
z_          =   -20*nm;
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
plot(z/nm,real(V_e),'k','LineWidth',1)
plot(z/nm,imag(V_e),'--k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$V^{e}_{i}$','Interpret','Latex')
title('Voltage','Interpret','Latex')
legend('Re','Im','Interpreter','Latex')
xlim([zN z0]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ve' -pdf -transparent
%% Current e 
figure()
hold on
plot(z/nm,real(I_e),'k','LineWidth',1)
plot(z/nm,imag(I_e),'--k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$I^{e}_{i}$','Interpret','Latex')
title('Current','Interpret','Latex')
legend('Re','Im','Interpreter','Latex')
xlim([zN z0]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ie' -pdf -transparent
%% Voltage h
figure()
hold on
plot(z/nm,real(V_h),'k','LineWidth',1)
plot(z/nm,imag(V_h),'--k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$V^{h}_{i}$','Interpret','Latex')
title('Voltage','Interpret','Latex')
legend('Re','Im','Interpreter','Latex')
xlim([zN z0]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Vh' -pdf -transparent
%% Current h 
figure()
hold on
plot(z/nm,real(I_h),'k','LineWidth',1)
plot(z/nm,imag(I_h),'--k','LineWidth',1)
y_range     =   PlotLayers(nm);
hold off
xlabel('$z$ [nm]','Interpret','Latex')
ylabel('$I^{h}_{i}$','Interpret','Latex')
title('Current','Interpret','Latex')
legend('Re','Im','Interpreter','Latex')
xlim([zN z0]/nm)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Ih' -pdf -transparent
%%