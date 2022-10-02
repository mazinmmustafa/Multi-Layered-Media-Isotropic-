close all; clear; clc;
%%
j           =   sqrt(-1);
%% Create Box
[~,lambda0,~,~,~,~,~]=Configs();
xc          =   0;
yc          =   0;
zc          =   -0.3*lambda0;
Lx          =   0.3*lambda0;
Ly          =   0.3*lambda0;
Lz          =   0.1*lambda0;
epsn        =   6-j*0.5;
option      =   2; % (1: Coarse) (2: Standard) (3: Fine)
%%
R        	=   0.5*lambda0;
Mesh=CreateBox(xc,yc,zc,Lx,Ly,Lz,epsn,option);
[N,~]       =   size(Mesh);
fprintf('N\t=\t%d\n',N);
%% Plot Shape
figure()
hold on
[N,~]      	=   size(Mesh);
for i=1:N
    xn          =   Mesh(i,1);
 	yn          =   Mesh(i,2);
  	zn          =   Mesh(i,3);
    dn          =   Mesh(i,4);
    vert = [(xn-dn/2)/lambda0 (yn-dn/2)/lambda0 (zn-dn/2)/lambda0;(xn+dn/2)/lambda0 ...
        (yn-dn/2)/lambda0 (zn-dn/2)/lambda0;(xn+dn/2)/lambda0 (yn+dn/2)/lambda0 ...
        (zn-dn/2)/lambda0;(xn-dn/2)/lambda0 (yn+dn/2)/lambda0 (zn-dn/2)/lambda0;...
        (xn-dn/2)/lambda0 (yn-dn/2)/lambda0 (zn+dn/2)/lambda0;(xn+dn/2)/lambda0 ...
        (yn-dn/2)/lambda0 (zn+dn/2)/lambda0;(xn+dn/2)/lambda0 (yn+dn/2)/lambda0 ...
        (zn+dn/2)/lambda0;(xn-dn/2)/lambda0 (yn+dn/2)/lambda0 (zn+dn/2)/lambda0];
    fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    patch('Vertices',vert,'Faces',fac,'FaceVertexCData',white(6),'FaceColor','flat')
end
plot3([-R R]/lambda0,[0 0],[0 0],'k','LineWidth',1)
plot3([0 0],[-R R]/lambda0,[0 0],'k','LineWidth',1)
plot3([0 0],[0 0],[-R R]/lambda0,'k','LineWidth',1)
hold off
xlabel('$x/\lambda_0$','Interpret','Latex')
ylabel('$y/\lambda_0$','Interpret','Latex')
zlabel('$z/\lambda_0$','Interpret','Latex')
title('Shape','Interpret','Latex')
axis equal
view(-45,20)
axis([-R R -R R -R R]/lambda0)
%%
