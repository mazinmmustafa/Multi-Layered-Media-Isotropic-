close all; clear; clc;
%%
j           =   sqrt(-1);
Ns          =   400;
%% Definitions
[~,lambda0,~,~,~,~,Data]=Configs();
y           =   0*lambda0;
z           =   0.5*lambda0;
x_min    	=   -2*lambda0;
x_max    	=   2*lambda0;
%% Corrections
x           =   linspace(x_min,x_max,Ns);
y           =   y*ones(size(x));
%% Load Data
ReJ         =   load(strcat(pwd,'\DDA\Data\ReJ.dat'));
ImJ         =   load(strcat(pwd,'\DDA\Data\ImJ.dat'));
J           =   ReJ+j*ImJ;
%% Load Mesh
ReMesh     	=   load(strcat(pwd,'\DDA\Data\ReMesh.dat'));
ImMesh    	=   load(strcat(pwd,'\DDA\Data\ImMesh.dat'));
Mesh      	=   ReMesh+j*ImMesh;
%% Compute Scattered Fields
Ex          =   zeros(1,Ns);
Ey          =   zeros(1,Ns);
Ez          =   zeros(1,Ns);
[N,~]      	=   size(Mesh);
tic;
count       =   0;
for i=1:Ns
    for P=1:N
        xm          =   Mesh(P,1);
        ym          =   Mesh(P,2);
        zm          =   Mesh(P,3);
        dm          =   Mesh(P,4);
        dVm         =   dm^3;
        [phi,rho]   =   cart2pol(x(i)-xm,y(i)-ym);  
        Jx          =   J(3*(P-1)+1);
        Jy          =   J(3*(P-1)+2);
        Jz          =   J(3*(P-1)+3);
        Ex(1,i)     =   Ex(1,i)+dVm*GEJxx(rho,phi,z,zm)*Jx+dVm*GEJxy(rho,phi,z,zm)*Jy+dVm*GEJxz(rho,phi,z,zm)*Jz;
        Ey(1,i)     =   Ey(1,i)+dVm*GEJyx(rho,phi,z,zm)*Jx+dVm*GEJyy(rho,phi,z,zm)*Jy+dVm*GEJyz(rho,phi,z,zm)*Jz;
        Ez(1,i)     =   Ez(1,i)+dVm*GEJzx(rho,phi,z,zm)*Jx+dVm*GEJzy(rho,phi,z,zm)*Jy+dVm*GEJzz(rho,phi,z,zm)*Jz;
        clc;
        fprintf('Step:\t%0.0f/100\n',count*100/(N*Ns));
        count       =   count+1;
    end
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%%
figure()
hold on
plot(x/lambda0,abs(Ex),'k','LineWidth',1)
plot(x/lambda0,abs(Ey),'--k','LineWidth',1)
plot(x/lambda0,abs(Ez),'-.k','LineWidth',1)
hold off
xlabel('$x/\lambda_0$','Interpret','Latex')
ylabel('$|E|$ [V/m]','Interpret','Latex')
title('Scattered Field','Interpret','Latex')
legend('$E_{x}$','$E_{y}$','$E_{z}$','Interpreter','Latex','Location','NorthEast')
xlim([x_min x_max]/lambda0)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\NearField' -pdf -transparent
%% Load FEKO Results 
Data_FEKO 	=   load(strcat(pwd,'\DDA\Data\FEKOResults1.dat'));
z_FEKO      =   Data_FEKO(:,1);
Ex_FEKO   	=   Data_FEKO(:,2);
Ey_FEKO   	=   Data_FEKO(:,3);
Ez_FEKO   	=   Data_FEKO(:,4);
%%
figure()
hold on
plot(x/lambda0,abs(Ex),'k','LineWidth',1)
plot(z_FEKO/lambda0,abs(Ex_FEKO),'--k','LineWidth',1)
plot(x/lambda0,abs(Ey),'b','LineWidth',1)
plot(z_FEKO/lambda0,abs(Ey_FEKO),'--b','LineWidth',1)
plot(x/lambda0,abs(Ez),'r','LineWidth',1)
plot(z_FEKO/lambda0,abs(Ez_FEKO),'--r','LineWidth',1)
hold off
xlabel('$x/\lambda_0$','Interpret','Latex')
ylabel('$|E|$ [V/m]','Interpret','Latex')
title('Scattered Field','Interpret','Latex')
legend('$E_{x}$ VIE','$E_{x}$ FEKO','$E_{y}$ VIE','$E_{y}$ FEKO','$E_{z}$ VIE','$E_{z}$ FEKO','Interpreter','Latex','Location','NorthEast')
xlim([x_min x_max]/lambda0)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\FEKO' -pdf -transparent
%%
% figure()
% hold on
% plot(x/lambda0,real(Ex),'k','LineWidth',1)
% plot(x/lambda0,imag(Ex),'--k','LineWidth',1)
% hold off
% xlabel('$x/\lambda_0$','Interpret','Latex')
% ylabel('$E$ [V/m]','Interpret','Latex')
% title('Scattered Field','Interpret','Latex')
% legend('$\textrm{Re}[E]$','$\textrm{Im}[E]$','Interpreter','Latex','Location','NorthEast')
% xlim([x_min x_max]/lambda0)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\NearField' -pdf -transparent
%% Compute Far Fields
% Nsf         =   1000;
% phi         =   0;
% theta       =   linspace(0,2*pi,Nsf);
% E_theta   	=   0;
% E_phi     	=   0;
% [N,~]      	=   size(Mesh);
% tic;
% count       =   0;
% for P=1:N
%   	xm          =   Mesh(P,1);
%   	ym          =   Mesh(P,2);
%   	zm          =   Mesh(P,3);
%  	dm          =   Mesh(P,4);
%   	dVm         =   dm^3;
%   	Jx          =   J(3*(P-1)+1);
%   	Jy          =   J(3*(P-1)+2);
%   	Jz          =   J(3*(P-1)+3);
%   	[E_theta_,E_phi_]=EFarField(theta,phi,Jx,Jy,Jz,xm,ym,zm);
%     E_theta     =   E_theta+dVm*E_theta_;
%     E_phi       =   E_phi+dVm*E_phi_;
%    	clc;
%   	fprintf('Step:\t%0.0f/100\n',count*100/N);
%   	count       =   count+1;
% end
% toc; clc;
% fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Plot
% figure()
% DATA_min  	=   -20;
% DATA        =   abs(E_theta)./max(abs(E_theta)); 
% DATA        =   20*log10(DATA);
% DATA(DATA<DATA_min)=DATA_min;
% PolardB(theta,DATA,[DATA_min 0],5,'k',1)
% title('$20\log_{10}|E_{\theta}|$ vs $\theta$','Interpret','Latex')
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\FarFieldTheta' -pdf -transparent
%% 