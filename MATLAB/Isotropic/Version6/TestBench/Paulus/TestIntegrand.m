close all; clear; clc;
%%
Ns          =   4000;
[~,~,~,~,nm,~,~,~,~]=Units();
%% Definitions
rho         =   633*nm;
phi         =   45;
z           =   750*nm;
z_          =   750*nm;
R           =   10;
%% Corrections
phi         =   deg2rad(phi);
kp          =   linspace(0,R,Ns);
%%
figure()
plot(kp,abs(Integrand(kp,rho,z,z_,0)),'k','LineWidth',1)
xlabel('$k_{\rho}/k_0$','Interpret','Latex')
ylabel('$|\tilde{G}^{EJ}|$','Interpret','Latex')
title('Integrand','Interpret','Latex')
% ylim([0 20e9])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\Integrand' -pdf -transparent
%%
function[G]=Integrand(kp,rho,z,z_,v)
[factor,kp]=Detour(kp,rho,abs(z-z_));
[V_e,~,V_h,~]=TLGFr(kp,'i',z,z_,1); % Without Direct Field
% [V_e,~,V_h,~]=TLGF(kp,'i',z,z_,1); % With Direct Field
G           =   factor.*(V_e+V_h).*besselj(v,kp*rho).*kp;   
end
%%
function[factor,kp]=Detour(kp,rho,Distance)
%% Substitution
[~,~,k0,~,~,~,~]=Configs();
kp         	=   k0*kp;
%% Detour
j          	=   sqrt(-1);
a          	=   k0*(Findkmin+1);
if rho>Distance
    b          	=   min(k0,1/rho);
else
    b           =   k0;
end
% b           =   0*k0;
t          	=   kp;
x         	=   t;
y         	=   b*sin(pi*t/a).*(t<a)+0.*(t>=a);
kp        	=   x+j*y;
factor    	=   1+j*(pi*b/a)*cos(pi*t/a); 
end
%%