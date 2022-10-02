function[Eu,Ev,Ez]=ModalFields(kp,z,Sheet)
%% Load Data
[~,~,~,omega,~,~,Data]=Configs();
[~,~,eps0,~,~,~,~]=Constants();
%%
[V_e,I_e,V_h,~]=TLGFm(kp,z,Sheet);
epsn        =   Data(SelectLayer(z)+1,3);
Eu          =   V_e;
Ev          =   V_h;
Ez          =   -kp*I_e./(omega*eps0*epsn);
end
function[V_e,I_e,V_h,I_h]=TLGFm(kp,z,Sheet)
j           =   sqrt(-1);
%% Find Layer
n           =   SelectLayer(z);
m           =   1;
%% Load Data
[~,~,~,~,~,~,Data]=Configs();
%% Solve for n=1
[kzn,~,Zen,Zhn]=kpParameters(kp,m,Sheet);
zn          = 	Data(m+1,1);
[Gamma_e_U,Gamma_h_U]=Refl(kp,'U',m,Sheet);
Vep         =   1;
Vhp         =   1;
%% Case n==m
if n==m
    V_e         =   Vep.*(exp(-j*kzn*(z-zn))+Gamma_e_U.*exp(j*kzn*(z-zn)));
  	V_h         =   Vhp.*(exp(-j*kzn*(z-zn))+Gamma_h_U.*exp(j*kzn*(z-zn)));
    I_e         =   Vep.*(exp(-j*kzn*(z-zn))-Gamma_e_U.*exp(j*kzn*(z-zn)))./Zen;
    I_h         =   Vhp.*(exp(-j*kzn*(z-zn))-Gamma_h_U.*exp(j*kzn*(z-zn)))./Zhn;
end
%%
[~,Theta_n,~,~]=kpParameters(kp,m+1,Sheet);
[Gamma_e_D2,Gamma_h_D2]=Refl(kp,'D',m+1,Sheet);
Vn_e      	=   Vep.*(1+Gamma_e_U)./(1+Gamma_e_D2.*exp(-j*2*Theta_n));
Vn_h      	=   Vhp.*(1+Gamma_h_U)./(1+Gamma_h_D2.*exp(-j*2*Theta_n));
%% Case n>m
if n>m
    [kzn,Theta_n,Zen,Zhn]=kpParameters(kp,n,Sheet);
    zn          = 	Data(n+1,1);
    [Gamma_e_D,Gamma_h_D]=Refl(kp,'D',n,Sheet);
    [Vn_e,Vn_h]=tau(kp,Vn_e,Vn_h,'D',m,n,Sheet);
    V_e         =   Vn_e.*exp(-j*Theta_n).*(exp(j*kzn*(z-zn))+Gamma_e_D.*exp(-j*kzn*(z-zn)));  
    V_h         =   Vn_h.*exp(-j*Theta_n).*(exp(j*kzn*(z-zn))+Gamma_h_D.*exp(-j*kzn*(z-zn)));
    I_e         =   -Vn_e.*exp(-j*Theta_n).*(exp(j*kzn*(z-zn))-Gamma_e_D.*exp(-j*kzn*(z-zn)))./Zen;
    I_h         =   -Vn_h.*exp(-j*Theta_n).*(exp(j*kzn*(z-zn))-Gamma_h_D.*exp(-j*kzn*(z-zn)))./Zhn;
end
end
%%