close all; clear; clc;
%%
j           =   sqrt(-1);
Ns          =   200;
%% Definitions
[~,lambda0,~,omega,~,~,Data]=Configs();
% Incident Plane Wave
phi_i       =   180;
theta_i     =   0;
E_TM        =   0;
E_TE        =   1;
% Box Center
xc          =   0;
yc          =   0;
zc          =   -0.5*lambda0+0.05*lambda0/2;
% Box Dimensions
Lx          =   0.3*lambda0;
Ly          =   0.3*lambda0;
Lz          =   0.05*lambda0;
% Scatterer Dielectric Constant
eps_s       =   10-j*5;
% Mesh Size
option      =   1; % (1: Coarse) (2: Standard) (3: Fine)
% Observation Points
x_min    	=   -2*lambda0;
x_max    	=   2*lambda0;
y           =   0*lambda0;
z           =   [0 0.1]*lambda0;
%% Corrections
x           =   linspace(x_min,x_max,Ns);
y           =   y*ones(Ns);
%% Create Mesh
Mesh=CreateBox(xc,yc,zc,Lx,Ly,Lz,eps_s,option);
[N,~]       =   size(Mesh);
fprintf('N\t=\t%d\n',N);
str         =   input('Enter to continue...');
%% Solve DDA
[J]=DDAPlaneWave(phi_i,theta_i,E_TM,E_TE,Mesh);
%%
Ey1      	=   zeros(1,Ns);
Ey2      	=   zeros(1,Ns);
tic;
count       =   0;
for i=1:Ns
    [~,Ey_,~]=ComputeFields(x(i),y(i),z(1),J,Mesh);
    Ey1(1,i) 	=   Ey1(1,i)+Ey_;
    [~,Ey_,~]=ComputeFields(x(i),y(i),z(2),J,Mesh);
    Ey2(1,i) 	=   Ey2(1,i)+Ey_;
    clc;
    fprintf('Computing Scattered Fields...\n');
	fprintf('Step:\t%0.0f/100\n',count*100/Ns);
	count       =   count+1;
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%%
figure()
hold on
plot(x/lambda0,abs(Ey1*1000),'k','LineWidth',1)
plot(x/lambda0,abs(Ey2*1000),'--k','LineWidth',1)
hold off
xlabel('$x/\lambda_0$','Interpret','Latex')
ylabel('$|E_{y}|$ [mV/m]','Interpret','Latex')
title('Scattered Field','Interpret','Latex')
xlim([x_min x_max]/lambda0)
legend('$z=0\lambda_0$','$z=0.1\lambda_0$','Interpreter','Latex','Location','NorthEast')
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\ScatteredFields' -pdf -transparent
%%