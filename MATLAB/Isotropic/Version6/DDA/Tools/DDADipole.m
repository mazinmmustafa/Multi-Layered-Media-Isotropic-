function[J]=DDADipole(x_,y_,z_,theta0,phi0,J,Mesh)
%%
j           =   sqrt(-1);
%% Constants
[~,~,eps0,eta0,~,~,~]=Constants();
[~,~,k0,omega,~,~,Data]=Configs();
%% Corrections
theta0      =   deg2rad(theta0);
phi0        =   deg2rad(phi0);
Jx          =   J*sin(theta0)*cos(phi0);
Jy          =   J*sin(theta0)*sin(phi0);
Jz          =   J*cos(theta0);
%% Create Mesh
[N,~]      	=   size(Mesh);
%% Create Matrix
count       =   0;
Z           =   zeros(3*N,3*N);
tic;
for P=1:N
    for Q=1:N
        fprintf('Creating Matrix...\n');
        fprintf('Step:\t%0.0f/100\n',100*count/N^2);
        dn          =   Mesh(Q,4);
        zn          =   Mesh(Q,3);
        epsn        =	Data(SelectLayer(zn)+1,3);
        kn          =	Data(SelectLayer(zn)+1,4);
        dVn         =   dn^3;
        an          =   dn*(3/(4*pi))^(1/3);
        M           =   (2/(3*kn^2))*((1+j*kn*an)*exp(-j*kn*an)-1);
        if P==Q
            zm          =   Mesh(P,3);
            epsm        =	Mesh(P,5);
            depsm       =   epsm-epsn;   
            term        =   ((M*(3*(kn^2)*depsm)/(2*epsn+epsm))-1);
            term        =   term*(2*epsn+epsm)/(j*omega*3*eps0*epsn*depsm);
            Z(3*(P-1)+1,3*(Q-1)+1)  =   term+dVn*GEJrxx(0,0,zm,zn); 
            Z(3*(P-1)+1,3*(Q-1)+2)  =   0;
            Z(3*(P-1)+1,3*(Q-1)+3)  =   0;
            Z(3*(P-1)+2,3*(Q-1)+1)  =   0;
            Z(3*(P-1)+2,3*(Q-1)+2)  =   term+dVn*GEJryy(0,0,zm,zn);
            Z(3*(P-1)+2,3*(Q-1)+3)  =   0;
            Z(3*(P-1)+3,3*(Q-1)+1)  =   0;
            Z(3*(P-1)+3,3*(Q-1)+2)  =   0;
            Z(3*(P-1)+3,3*(Q-1)+3)  =   term+dVn*GEJrzz(0,0,zm,zn)-(eta0/(j*k0*epsn));    
        end
        if P~=Q
            xm          =   Mesh(P,1);
            ym          =   Mesh(P,2);
            zm          =   Mesh(P,3);
            xn          =   Mesh(Q,1);
            yn          =   Mesh(Q,2);
            zn          =   Mesh(Q,3);
            [phi_nm,rho_mn]   =   cart2pol(xn-xm,yn-ym);
            Z(3*(P-1)+1,3*(Q-1)+1)  =   dVn*GEJxx(rho_mn,phi_nm,zm,zn); 
            Z(3*(P-1)+1,3*(Q-1)+2)  =   dVn*GEJxy(rho_mn,phi_nm,zm,zn); 
            Z(3*(P-1)+1,3*(Q-1)+3)  =   dVn*GEJxz(rho_mn,phi_nm,zm,zn); 
            Z(3*(P-1)+2,3*(Q-1)+1)  =   dVn*GEJyx(rho_mn,phi_nm,zm,zn); 
            Z(3*(P-1)+2,3*(Q-1)+2)  =   dVn*GEJyy(rho_mn,phi_nm,zm,zn); 
            Z(3*(P-1)+2,3*(Q-1)+3)  =   dVn*GEJyz(rho_mn,phi_nm,zm,zn); 
            Z(3*(P-1)+3,3*(Q-1)+1)  =   dVn*GEJzx(rho_mn,phi_nm,zm,zn); 
            Z(3*(P-1)+3,3*(Q-1)+2)  =   dVn*GEJzy(rho_mn,phi_nm,zm,zn); 
            Z(3*(P-1)+3,3*(Q-1)+3)  =   dVn*GEJzz(rho_mn,phi_nm,zm,zn); 
        end
        count       =   count+1;
        clc;
    end
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Incident Field
Ei          =   zeros(3*N,1);
clc;
count       =   0;
for P=1:N
    fprintf('Computing Incident Fields...\n');
    fprintf('Step:\t%0.0f/100\n',100*count/N);
    xm          =   Mesh(P,1);
  	ym          =   Mesh(P,2);
   	zm          =   Mesh(P,3);
    [phi_m,rho_m]   =   cart2pol(xm-x_,ym-y_);
  	Ei(3*(P-1)+1,1)	=   GEJxx(rho_m,phi_m,zm,z_)*Jx+GEJxy(rho_m,phi_m,zm,z_)*Jy+GEJxz(rho_m,phi_m,zm,z_)*Jz;
  	Ei(3*(P-1)+2,1)	=   GEJyx(rho_m,phi_m,zm,z_)*Jx+GEJyy(rho_m,phi_m,zm,z_)*Jy+GEJyz(rho_m,phi_m,zm,z_)*Jz;
	Ei(3*(P-1)+3,1)	=   GEJzx(rho_m,phi_m,zm,z_)*Jx+GEJzy(rho_m,phi_m,zm,z_)*Jy+GEJzz(rho_m,phi_m,zm,z_)*Jz;
    count       =   count+1;
    clc;
end
%% Compute J
J           =   -Z\Ei;
end
%%