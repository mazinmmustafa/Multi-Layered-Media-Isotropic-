close all; clear; clc;
%%
addpath('colormap')
%%
um      =   1E-6;
%%
Axis    =   load('Axis.dat');
Ex      =   load('DataEx.dat');
Ey      =   load('DataEy.dat');
Ez      =   load('DataEz.dat');
%%
[x,z]  	=   meshgrid(Axis(:,1),Axis(:,2));
%%
figure()
pcolor(x/um,z/um,20*log10(abs(Ex)'))
hold on
plot([-2 +2],[0 0],'-k','LineWidth',1)
plot([-2 +2],[-0.05 -0.05],'-k','LineWidth',1)
hold off
colormap jet
shading flat
axis equal
colorbar
xlabel('$x$ [$\mu$m]','Interpret','Latex','FontSize',15)
ylabel('$z$ [$\mu$m]','Interpret','Latex','FontSize',15)
title('$20\log_{10}|\mathcal{R}e[E_{x}]|$','Interpret','Latex','FontSize',15)
set(gca,'TickLabel','Latex','FontSize',15)
set(colorbar,'TickLabelInterpreter','Latex','FontSize',15)
axis([-1 +1 -1 +1]*2)
caxis([220 340])
%%
% exportgraphics(gcf,'FigureEx.png','Resolution',200)
%%
figure()
pcolor(x/um,z/um,20*log10(abs(Ez)'))
hold on
plot([-2 +2],[0 0],'-k','LineWidth',1)
plot([-2 +2],[-0.05 -0.05],'-k','LineWidth',1)
hold off
colormap(twilight_shifted)
% colormap(twilight)
shading flat
axis equal
colorbar
xlabel('$x$ [$\mu$m]','Interpret','Latex','FontSize',15)
ylabel('$z$ [$\mu$m]','Interpret','Latex','FontSize',15)
% title('$20\log_{10}|\mathcal{R}e[E_{z}]|$','Interpret','Latex','FontSize',15)
set(gca,'TickLabel','Latex','FontSize',15)
set(colorbar,'TickLabelInterpreter','Latex','FontSize',15)
axis([-1 +1 -1 +1]*2)
caxis([220 340])
axis off
colorbar off
%%
exportgraphics(gcf,'FigureEz.png','Resolution',200)
%%

