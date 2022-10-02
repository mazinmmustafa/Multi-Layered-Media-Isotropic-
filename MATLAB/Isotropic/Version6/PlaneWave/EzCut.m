function[E_rho,E_phi,E_z]=EzCut(rho,phi,z,theta_i,phi_i,Ei_theta,Ei_phi)
%% Load Data
[~,~,~,~,~,~,Data]=Configs();
z0          =   Data(1,1);
eps1     	=   Data(2,3);
k1          =   Data(2,4);
kpi         =   k1*sin(theta_i);
n           =   SelectLayer(z);
[~,~,Ze1,~]=kpParameters(kpi,1,1);
epsn        =   Data(n+1,3);   
%%
j           =   sqrt(-1);
factor      =   -2*exp(j*kpi*rho.*cos(phi_i-phi)).*exp(j*k1*cos(theta_i)*z0);
[V_e,I_e,V_h,~]=TLGF(kpi,'v',z,z0,1);
E_rho     	= 	Ei_theta*factor.*V_e*cos(theta_i);
E_phi       =   Ei_phi*factor.*V_h;       
E_z         =   -Ei_theta*factor.*I_e.*Ze1*(eps1/epsn)*sin(theta_i);
end
%%