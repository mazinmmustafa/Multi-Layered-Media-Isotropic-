function[]=PlotLayersH(x_min,x_max,unit)
%% Load Data
[N,~,~,~,~,~,Data]=Configs();
%% Plot Lines
for n=1:N-1  
    plot([x_min x_max]/unit,[Data(n+1,1) Data(n+1,1)]/unit,'k','LineWidth',0.5)
end
end
%%