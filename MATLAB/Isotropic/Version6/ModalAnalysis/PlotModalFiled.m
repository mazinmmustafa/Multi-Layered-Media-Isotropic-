close all; clear; clc;
%%
[~,~,~,um,~,~,~,~,~]=Units();
j           =   sqrt(-1);
Ns          =   1000;
%%
Sheet       =   1;
kp          =   1.88224222918665-j*0.00063471402154; pol='e';
% kp          =   3.50344333295000+j*0.00710300097870; pol='h';
% kp          =   3.33728685820780-j*0.00022949110400; pol='h';
% kp          =   3.25168520698340-j*0.00053051477990; pol='h';
% kp          =   3.10425142141457+j*0.00133798633975; pol='h';
% kp          =   2.62813932045903+j*0.00154864433115; pol='h';
% kp          =   1.76819096041243+j*0.00135321718386; pol='h';
% kp          =   1.07426202652578+j*0.00245789147357; pol='h';
% kp          =   3.49668379589130+j*0.00654398171100; pol='e';
% kp          =   3.33069711910720+j*0.00003518642230; pol='e';
% kp          =   3.22433799874650-j*0.00017448261260; pol='e';
% kp          =   2.79439777568252+j*0.00070878520448; pol='e';
% kp          =   2.46292446281425+j*0.00117932006477; pol='e';
% kp          =   2.00514007332263+j*0.00160292202929; pol='e';
% kp          =   1.35099878658162+j*0.00231404951497; pol='e';
% kp          =   1.00143843982593+j*0.00004669412354; pol='e';
%% 
[N,~,k0,~,~,~,Data]=Configs();  
z0        	=   Data(1,1);
zN        	=   Data(N+1,1);
D           =   z0-zN;
% z           =   linspace(zN-D/2,z0+D/2,Ns);
z           =   linspace(-200*um,150*um,Ns); % Manual limits
kp          =   kp*k0;
%%
Eu          =   zeros(1,Ns);
Ev          =   zeros(1,Ns);
Ez          =   zeros(1,Ns);
for i=1:Ns
    [Eu(1,i),Ev(1,i),Ez(1,i)]=ModalFields(kp,z(i),Sheet);
end
%%
if pol=='e'
    E           =   sqrt(abs(Eu).^2+abs(Ez).^2);
end
if pol=='h'
    E           =   abs(Ev);
end
%% Normalized
E           =   E./max(E);
%% E Field
figure()
hold on
plot(z/um,E,'k','LineWidth',1)
y_range     =   PlotLayers(um);
hold off
xlabel('$z$ [$\mu$m]','Interpret','Latex')
ylabel('$|E|$','Interpret','Latex')
title('Field Profile','Interpret','Latex')
xlim([min(z) max(z)]/um)
ylim(y_range)
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\EProfile' -pdf -transparent
%%
