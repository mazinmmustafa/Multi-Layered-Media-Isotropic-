close all; clear; clc;
%%
j           =   sqrt(-1);
Nc          =   50;
Ns          =   200;
%% Definitions
pol         =   'e';
R           =   5;
Sheet       =   1;
%%
[x,y]       =   meshgrid(linspace(-R,R,Ns),linspace(-R,R,Ns));
kp          =   x+j*y;
%%
DF          =   TMatrix(kp,pol,Sheet);
%% Plot 
figure()
pcolor(x,y,log10(abs(DF)))
hold on
contour(x,y,log10(abs(DF)),Nc,'k')
plot([-R R],[0 0],'k','LineWidth',1)
plot([0 0],[-R R],'k','LineWidth',1)
hold off
colormap jet
shading flat
colorbar 
xlabel('$k^{\prime}_{\rho}/k_{0}$','Interpret','Latex')
ylabel('$k^{\prime\prime}_{\rho}/k_{0}$','Interpret','Latex')
title('$\log_{10}|\textrm{DF}|$','Interpret','Latex')
axis([-R R -R R])
% caxis([2 14])
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\DF' -pdf -transparent
%% GNU plot
% gnuplot_pcolor(x,y,log10(abs(DF))),'Data_Matrix.dat');
%%

