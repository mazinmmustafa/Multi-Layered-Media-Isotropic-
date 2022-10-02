close all; clear; clc;
%%
j           =   sqrt(-1);
%% Constants
[~,~,eps0,~,~,~,~]=Constants();
[~,lambda0,~,omega,~,~,Data]=Configs();
%% Definitions
% Plane Wave
theta_i     =   0;
phi_i       =   0;
E_TM        =   0;
E_TE        =   1;
% Dipole
Jx          =   1;
Jy          =   0;
Jz          =   0;
x_          =   0.4*lambda0;
y_          =   0.6*lambda0;
z_          =   0.5*lambda0;
% Box
xc          =   0;
yc          =   0;
zc          =   -0.35*lambda0;
Lx          =   0.3*lambda0;
Ly          =   0.3*lambda0;
Lz          =   0.05*lambda0;
eps_s       =   10-j*5;
option      =   1; % (1: Coarse) (2: Standard) (3: Fine)
%% Corrections
theta_i     =   deg2rad(theta_i);
phi_i       =   deg2rad(phi_i);
%% Create Mesh
Mesh=CreateBox(xc,yc,zc,Lx,Ly,Lz,eps_s,option);
[N,~]      	=   size(Mesh);
%% Save Mesh
DATA = real(Mesh); save(strcat(pwd,'\DDA\Data\','ReMesh.dat'),'DATA','-ascii');
DATA = imag(Mesh); save(strcat(pwd,'\DDA\Data\','ImMesh.dat'),'DATA','-ascii');
%% Create Matrix
count       =   0;
Z           =   zeros(3*N,3*N);
tic;
for P=1:N
    for Q=1:N
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
            Z(3*(P-1)+1,3*(Q-1)+2)  =   dVn*GEJrxy(0,0,zm,zn);
            Z(3*(P-1)+1,3*(Q-1)+3)  =   dVn*GEJrxz(0,0,zm,zn);
            Z(3*(P-1)+2,3*(Q-1)+1)  =   dVn*GEJryx(0,0,zm,zn);
            Z(3*(P-1)+2,3*(Q-1)+2)  =   term+dVn*GEJryy(0,0,zm,zn);
            Z(3*(P-1)+2,3*(Q-1)+3)  =   dVn*GEJryz(0,0,zm,zn);
            Z(3*(P-1)+3,3*(Q-1)+1)  =   dVn*GEJrzx(0,0,zm,zn);
            Z(3*(P-1)+3,3*(Q-1)+2)  =   dVn*GEJrzy(0,0,zm,zn);
            Z(3*(P-1)+3,3*(Q-1)+3)  =   term+dVn*GEJrzz(0,0,zm,zn);     
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
%% Plane Wave
% for P=1:N
%     xm          =   Mesh(P,1);
%   	ym          =   Mesh(P,2);
%    	zm          =   Mesh(P,3);
%     [Ex,Ey,Ez]=EiPlaneWave(xm,ym,zm,theta_i,phi_i,E_TM,E_TE);
%     Ei(3*(P-1)+1,1) =   Ex;
%     Ei(3*(P-1)+2,1) =   Ey;
%     Ei(3*(P-1)+3,1) =   Ez;
% end
%% Dipole
count       =   0;
for P=1:N
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
%% Save Data
DATA = real(J); save(strcat(pwd,'\DDA\Data\','ReJ.dat'),'DATA','-ascii');
DATA = imag(J); save(strcat(pwd,'\DDA\Data\','ImJ.dat'),'DATA','-ascii');
%% 