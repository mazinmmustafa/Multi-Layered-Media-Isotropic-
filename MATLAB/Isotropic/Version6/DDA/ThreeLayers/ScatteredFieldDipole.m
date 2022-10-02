close all; clear; clc;
%%
j           =   sqrt(-1);
Ns          =   50;
%% Definitions
[~,lambda0,~,omega,~,~,Data]=Configs();
% Incident Field
x_          =   0.4*lambda0;
y_          =   0.6*lambda0;
z_          =   0.5*lambda0;
theta0      =   60;
phi0        =   135;
J           =   1;
% Box Center
xc          =   0;
yc          =   0;
zc          =   -0.3*lambda0;
% Box Dimensions
Lx          =   0.3*lambda0;
Ly          =   0.3*lambda0;
Lz          =   0.1*lambda0;
% Scatterer Dielectric Constant
eps_s       =   6-j*0.5;
% Mesh Size
option      =   2; % (1: Coarse) (2: Standard) (3: Fine)
% Observation Points
x_min    	=   -2*lambda0;
x_max    	=   2*lambda0;
y           =   0*lambda0;
z           =   0.5*lambda0;
%% Corrections
x           =   linspace(x_min,x_max,Ns);
y           =   y*ones(Ns);
%% Create Mesh
Mesh=CreateBox(xc,yc,zc,Lx,Ly,Lz,eps_s,option);
[N,~]       =   size(Mesh);
fprintf('N\t=\t%d\n',N);
str         =   input('Enter to continue...');
%% Solve DDA
tic;
[J]=DDADipole(x_,y_,z_,theta0,phi0,J,Mesh);
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Compute Scattered Fields
Ex          =   zeros(1,Ns);
Ey          =   zeros(1,Ns);
Ez          =   zeros(1,Ns);
[N,~]      	=   size(Mesh);
tic;
count       =   0;
for i=1:Ns
    [Ex_,Ey_,Ez_]=ComputeFields(x(i),y(i),z,J,Mesh);
    Ex(1,i)     =   Ex(1,i)+Ex_;
    Ey(1,i)     =   Ey(1,i)+Ey_;
    Ez(1,i)     =   Ez(1,i)+Ez_;
    clc;
	fprintf('Computing Scattered Fields\n');
 	fprintf('Step:\t%0.0f/100\n',count*100/Ns);
	count       =   count+1;
end
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%% Load FEKO Results 
Data_FEKO 	=   load(strcat(pwd,'\DDA\Data\FEKOThreeLayers0_5lambda0.dat'));
z_FEKO      =   Data_FEKO(:,1);
Ex_FEKO   	=   Data_FEKO(:,2);
Ey_FEKO   	=   Data_FEKO(:,3);
Ez_FEKO   	=   Data_FEKO(:,4);
%%
figure()
hold on
plot(z_FEKO/lambda0,abs(Ex_FEKO),'k','LineWidth',1)
plot(x/lambda0,abs(Ex),'ok','LineWidth',1)
plot(z_FEKO/lambda0,abs(Ey_FEKO),'b','LineWidth',1)
plot(x/lambda0,abs(Ey),'ob','LineWidth',1)
plot(z_FEKO/lambda0,abs(Ez_FEKO),'r','LineWidth',1)
plot(x/lambda0,abs(Ez),'or','LineWidth',1)
hold off
xlabel('$x/\lambda_0$','Interpret','Latex')
ylabel('$|E|$ [V/m]','Interpret','Latex')
title('Scattered Field','Interpret','Latex')
legend('$E_{x}$ FEKO','$E_{x}$ VIE','$E_{y}$ FEKO','$E_{y}$ VIE','$E_{z}$ FEKO','$E_{z}$ VIE','Interpreter','Latex','Location','NorthEast')
xlim([x_min x_max]/lambda0)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\FEKO' -pdf -transparent
%%