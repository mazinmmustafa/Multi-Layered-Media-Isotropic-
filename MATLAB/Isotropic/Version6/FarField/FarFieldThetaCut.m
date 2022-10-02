function[E_theta,E_phi]=FarFieldThetaCut(theta,phi,theta0,phi0,x_,y_,z_,J)
Ns          =   length(theta);
j           =   sqrt(-1);
%% Corrections
phi         =   deg2rad(phi);  
theta0      =   deg2rad(theta0);
phi0        =   deg2rad(phi0);
Jx          =   J*sin(theta0)*cos(phi0);
Jy          =   J*sin(theta0)*sin(phi0);
Jz          =   J*cos(theta0);
phi         =   phi*ones(1,Ns);
for i=1:Ns
    if theta(1,i)>pi
        theta(1,i)  =   2*pi-theta(1,i);
        phi(1,i) 	=   phi(1,i)+pi;
    end
end
%% Load Data
[N,~,~,~,~,~,Data]=Configs();
z0          =   Data(1,1);
eps1     	=   Data(2,3);
k1          =   Data(2,4);
eta1     	=   Data(2,5);
zN          =   Data(N+1,1);
epsN     	=   Data(N+1,3);
kN          =   Data(N+1,4);
etaN     	=   Data(N+1,5);
m           =   SelectLayer(z_);
eps_        =   Data(m+1,3);   
%% Up
kp1      	=   k1*sin(theta);
factor      =   -j*k1*2*exp(j*k1*cos(theta)*z0).*exp(j*k1*sin(theta).*(x_*cos(phi)+y_*sin(phi)));      
[Vei,~,Vhi,~]=TLGF(kp1,'i',z0,z_,1);
[Vev,~,~,~]=TLGF(kp1,'v',z0,z_,1);
E_theta     =   Vei.*(cos(phi)*Jx+sin(phi)*Jy)-Vev.*(eta1*eps1/eps_).*sin(theta)*Jz;
E_phi       =   Vhi.*cos(theta).*(-sin(phi)*Jx+cos(phi)*Jy);
E_theta_U  	=   E_theta.*factor;
E_phi_U    	=   E_phi.*factor;
%% Down
kpN      	=   kN*sin(theta);
factor      =   -j*kN*2*exp(j*kN*cos(theta)*zN).*exp(j*kN*sin(theta).*(x_*cos(phi)+y_*sin(phi)));      
[Vei,~,Vhi,~]=TLGF(kpN,'i',zN,z_,1);
[Vev,~,~,~]=TLGF(kpN,'v',zN,z_,1);
E_theta     =   Vei.*(cos(phi)*Jx+sin(phi)*Jy)-Vev.*(etaN*epsN/eps_).*sin(theta)*Jz;
E_phi       =   Vhi.*cos(theta).*(-sin(phi)*Jx+cos(phi)*Jy);
E_theta_D  	=   E_theta.*factor;
E_phi_D    	=   E_phi.*factor;
%% Total
E_theta     =   E_theta_U.*(theta<=pi/2)+E_theta_D.*(theta>pi/2); 
E_phi       =   E_phi_U.*(theta<=pi/2)+E_phi_D.*(theta>pi/2); 
%% Remove Small Values
E_theta(abs(E_theta)<1e-2)=0;
E_phi(abs(E_phi)<1e-2)=0;
end