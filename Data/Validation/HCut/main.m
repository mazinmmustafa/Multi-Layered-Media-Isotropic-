close all; clear; clc;
%%
Data 	= 	load('Data1.dat');
EJx     =   load('EJx.dat');
EJy     =   load('EJy.dat');
EJz     =   load('EJz.dat');
HJx     =   load('HJx.dat');
HJy     =   load('HJy.dat');
HJz     =   load('HJz.dat');
%%
figure()
hold on
plot(Data(:,1),Data(:,2),'-k','LineWidth',1)
plot(EJx(:,1),EJx(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|E_{x}|$ [V/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 1400])
%%
exportgraphics(gcf,'FigureEJx.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,3),'-k','LineWidth',1)
plot(EJy(:,1),EJy(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|E_{y}|$ [V/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 1400])
%%
exportgraphics(gcf,'FigureEJy.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,4),'-k','LineWidth',1)
plot(EJz(:,1),EJz(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|E_{z}|$ [V/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 3000])
%%
exportgraphics(gcf,'FigureEJz.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,5),'-k','LineWidth',1)
plot(HJx(:,1),HJx(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|H_{x}|$ [A/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 4])
%%
exportgraphics(gcf,'FigureHJx.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,6),'-k','LineWidth',1)
plot(HJy(:,1),HJy(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|H_{y}|$ [A/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 6])
%%
exportgraphics(gcf,'FigureHJy.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,7),'-k','LineWidth',1)
plot(HJz(:,1),HJz(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|H_{z}|$ [A/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 2])
%%
exportgraphics(gcf,'FigureHJz.pdf','ContentType','vector')
%%
Data 	= 	load('Data2.dat');
EMx     =   load('EMx.dat');
EMy     =   load('EMy.dat');
EMz     =   load('EMz.dat');
HMx     =   load('HMx.dat');
HMy     =   load('HMy.dat');
HMz     =   load('HMz.dat');
%%
figure()
hold on
plot(Data(:,1),Data(:,2),'-k','LineWidth',1)
plot(EMx(:,1),EMx(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|E_{x}|$ [V/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 6])
%%
exportgraphics(gcf,'FigureEMx.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,3),'-k','LineWidth',1)
plot(EMy(:,1),EMy(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|E_{y}|$ [V/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 6])
%%
exportgraphics(gcf,'FigureEMy.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,4),'-k','LineWidth',1)
plot(EMz(:,1),EMz(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|E_{z}|$ [V/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
xlim([-3 +3])
ylim([0 2])
%%
exportgraphics(gcf,'FigureEMz.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,5),'-k','LineWidth',1)
plot(HMx(:,1),HMx(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|H_{x}|$ [A/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
ax      =	gca;
ax.YAxis.Exponent = -3;
xlim([-3 +3])
ylim([0 8E-3])
%%
exportgraphics(gcf,'FigureHMx.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,6),'-k','LineWidth',1)
plot(HMy(:,1),HMy(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|H_{y}|$ [A/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
ax      =	gca;
ax.YAxis.Exponent = -3;
xlim([-3 +3])
ylim([0 5E-3])
%%
exportgraphics(gcf,'FigureHMy.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data(:,1),Data(:,7),'-k','LineWidth',1)
plot(HMz(:,1),HMz(:,2),'ok')
hold off
xlabel('$x$ [m]','Interpret','Latex','FontSize',14)
ylabel('$|H_{z}|$ [A/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('SI','FEKO','Interpreter','Latex','FontSize',14)
ax      =	gca;
ax.YAxis.Exponent = -3;
xlim([-3 +3])
ylim([0 12E-3])
%%
exportgraphics(gcf,'FigureHMz.pdf','ContentType','vector')
%%