close all; clear; clc;
%%
[N,lambda0,k0,omega,Gamma_U,Gamma_D,Data]=Configs();
%%
% kp      =   0.7*k0;
% [kz,Theta,Ze,Zh]=kpParameters(kp,3,1);
% fprintf('%21.14E, %21.14E\n',real(kz),imag(kz));
% fprintf('%21.14E, %21.14E\n',real(Theta),imag(Theta));
% fprintf('%21.14E, %21.14E\n',real(Ze),imag(Ze));
% fprintf('%21.14E, %21.14E\n',real(Zh),imag(Zh));
%%
x       =   lambda0;
y       =   lambda0/2;
z       =   -700E-9;
x_      =   0;
y_      =   0;
z_    	=   +300E-9;
rho     =   sqrt((x-x_)^2+(y-y_)^2);
phi     =   atan2(y-y_,x-x_);
GEJxx(rho,phi,z,z_)
GEJxy(rho,phi,z,z_)
GEJxz(rho,phi,z,z_)
GEJyx(rho,phi,z,z_)
GEJyy(rho,phi,z,z_)
GEJyz(rho,phi,z,z_)
GEJzx(rho,phi,z,z_)
GEJzy(rho,phi,z,z_)
GEJzz(rho,phi,z,z_)