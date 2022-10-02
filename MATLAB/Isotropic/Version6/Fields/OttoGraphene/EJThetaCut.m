close all; clear; clc;
%%
Ns          =   1000;
[~,~,~,um,~,~,~,~,~]=Units();
%% Definitions
[~,lambda0,~,eta0,~,~,~]=Configs();
theta0      =   0;
phi0        =   0;
J           =   1;
x_          =   0*um;
y_          =   0*um;
z_          =   5*um;
R           =   10*lambda0;
phi         =   0;
%% Corrections
phi_save    =   phi;
phi         =   deg2rad(phi);
phi         =   phi*ones(1,Ns);
theta       =   linspace(0,2*pi,Ns);
for i=1:Ns
    if theta(1,i)>pi
        theta(1,i)      =   2*pi-theta(1,i);
        phi(1,i)        =   phi(1,i)+pi;
    end
end
x           =   R*sin(theta).*cos(phi);
y           =   R*sin(theta).*sin(phi);
z           =   R*cos(theta);
[phi,rho]   =   cart2pol(x-x_,y-y_);
%% Dipole Components
[Jx,Jy,Jz]=DipoleJ(theta0,phi0);
%% Compute Fields
Ex          =   zeros(1,Ns);
Ey          =   zeros(1,Ns);
Ez          =   zeros(1,Ns);
tic;
count       =   0;
for i=1:Ns
	Ex(1,i)     =   GEJxx(rho(1,i),phi(1,i),z(1,i),z_)*Jx+GEJxy(rho(1,i),phi(1,i),z(1,i),z_)*Jy+GEJxz(rho(1,i),phi(1,i),z(1,i),z_)*Jz;
  	Ey(1,i)     =   GEJyx(rho(1,i),phi(1,i),z(1,i),z_)*Jx+GEJyy(rho(1,i),phi(1,i),z(1,i),z_)*Jy+GEJyz(rho(1,i),phi(1,i),z(1,i),z_)*Jz;
  	Ez(1,i)     =   GEJzx(rho(1,i),phi(1,i),z(1,i),z_)*Jx+GEJzy(rho(1,i),phi(1,i),z(1,i),z_)*Jy+GEJzz(rho(1,i),phi(1,i),z(1,i),z_)*Jz;
  	count       =   count+1;
  	clc;
 	fprintf('Step:\t%0.0f/100\n',count*100/Ns);
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Exact Fields
E_theat     =   Ex.*cos(theta).*cos(phi)+Ey.*cos(theta).*sin(phi)-Ez.*sin(theta);
E_phi       =   -Ex.*sin(phi)+Ey.*cos(phi);
%% Far Fields
[E_theta_far,E_phi_far]=FarFieldThetaCut(theta,phi_save,theta0,phi0,x_,y_,z_,J);
%% Plot E_theta
figure()
Data_min    =   -40;
Data_max    =   -8;
theta       =   linspace(0,2*pi,Ns);
Data        =   R*abs(E_theat)/sqrt(2*eta0); 
Data        =   20*log10(Data);
Data(Data<Data_min)=Data_min;
PolardB(theta,Data,[Data_min Data_max],5,'k',1)
hold on
Data        =   abs(E_theta_far)/(4*pi*sqrt(2*eta0)); 
Data        =   20*log10(Data);
Data(Data<Data_min)=Data_min;
PolardB(theta,Data,[Data_min Data_max],5,'--k',1)
hold off
title('$20\log_{10}|E_{\theta}|$ vs $\theta$','Interpret','Latex')
legend('Exact','Far Field','Interpreter','Latex','Location','NorthEast')
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\FarFieldTheta' -pdf -transparent
%% Plot E_phi
% figure()
% Data_min    =   -40;
% Data_max    =   0;
% theta       =   linspace(0,2*pi,Ns);
% Data        =   R*abs(E_phi)/sqrt(2*eta0); 
% Data        =   20*log10(Data);
% Data(Data<Data_min)=Data_min;
% PolardB(theta,Data,[Data_min Data_max],5,'k',1)
% hold on
% Data        =   abs(E_phi_far)/(4*pi*sqrt(2*eta0)); 
% Data        =   20*log10(Data);
% Data(Data<Data_min)=Data_min;
% PolardB(theta,Data,[Data_min Data_max],5,'--k',1)
% hold off
% title('$20\log_{10}|E_{\phi}|$ vs $\theta$','Interpret','Latex')
% legend('Exact','Far Field','Interpreter','Latex','Location','NorthEast')
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\FarFieldPhi' -pdf -transparent
%%
function[Jx,Jy,Jz]=DipoleJ(theta0,phi0)
theta0      =   deg2rad(theta0);
phi0        =   deg2rad(phi0);
Jx          =   sin(theta0)*cos(phi0);
Jy          =   sin(theta0)*sin(phi0);
Jz          =   cos(theta0);
end
%%










%%