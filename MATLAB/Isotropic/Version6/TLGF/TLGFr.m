function[V_e,I_e,V_h,I_h]=TLGFr(kp,s,z,z_,Sheet)
j           =   sqrt(-1);
%% Find Layers
m           =   SelectLayer(z_);
n           =   SelectLayer(z);
%% Load Data
[~,~,~,~,~,~,Data]=Configs();
%% Define Sources
if s=='v'
  	v           =   1;
   	i           =   0;
end
if s=='i'
 	v           =   0;
  	i           =   1;
end   
%% Solve for n=m
[kzn,Theta_n,Zen,Zhn]=kpParameters(kp,m,Sheet);
zn          = 	Data(m+1,1);
zm          = 	Data(m,1);
[Gamma_e_D,Gamma_h_D]=Refl(kp,'D',m,Sheet);
[Gamma_e_U,Gamma_h_U]=Refl(kp,'U',m,Sheet);
De          =   1-Gamma_e_D.*Gamma_e_U.*exp(-j*2*Theta_n);
Dh          =   1-Gamma_h_D.*Gamma_h_U.*exp(-j*2*Theta_n);
Vep         =   ((1-Gamma_e_D.*exp(-j*2*kzn*(z_-zn)))*v+(1+Gamma_e_D.*exp(-j*2*kzn*(z_-zn))).*Zen*i)./(2*De);
Ven         =   (-(1-Gamma_e_U.*exp(-j*2*kzn*(zm-z_)))*v+(1+Gamma_e_U.*exp(-j*2*kzn*(zm-z_))).*Zen*i)./(2*De);
Vhp         =   ((1-Gamma_h_D.*exp(-j*2*kzn*(z_-zn)))*v+(1+Gamma_h_D.*exp(-j*2*kzn*(z_-zn))).*Zhn*i)./(2*Dh);
Vhn         =   (-(1-Gamma_h_U.*exp(-j*2*kzn*(zm-z_)))*v+(1+Gamma_h_U.*exp(-j*2*kzn*(zm-z_))).*Zhn*i)./(2*Dh);
Gamma_e_U_  =   Gamma_e_U.*exp(-j*2*kzn*(zm-z_));
Gamma_e_D_  =   Gamma_e_D.*exp(-j*2*kzn*(z_-zn));
Gamma_h_U_  =   Gamma_h_U.*exp(-j*2*kzn*(zm-z_));
Gamma_h_D_  =   Gamma_h_D.*exp(-j*2*kzn*(z_-zn));
%% Case n=m
if n==m
    Gamma_e_D   =   0;
    Gamma_h_D   =   0;
    Gamma_e_U   =   0;
    Gamma_h_U   =   0;
    De0       	=   1-Gamma_e_D.*Gamma_e_U.*exp(-j*2*Theta_n);
    Dh0      	=   1-Gamma_h_D.*Gamma_h_U.*exp(-j*2*Theta_n);
    Vep0       	=   ((1-Gamma_e_D.*exp(-j*2*kzn*(z_-zn)))*v+(1+Gamma_e_D.*exp(-j*2*kzn*(z_-zn))).*Zen*i)./(2*De0);
    Ven0      	=   (-(1-Gamma_e_U.*exp(-j*2*kzn*(zm-z_)))*v+(1+Gamma_e_U.*exp(-j*2*kzn*(zm-z_))).*Zen*i)./(2*De0);
    Vhp0       	=   ((1-Gamma_h_D.*exp(-j*2*kzn*(z_-zn)))*v+(1+Gamma_h_D.*exp(-j*2*kzn*(z_-zn))).*Zhn*i)./(2*Dh0);
    Vhn0       	=   (-(1-Gamma_h_U.*exp(-j*2*kzn*(zm-z_)))*v+(1+Gamma_h_U.*exp(-j*2*kzn*(zm-z_))).*Zhn*i)./(2*Dh0);
    if z>=z_
        V_e         =   ((Vep-Vep0).*exp(-j*kzn*(z-z_))+Vep.*Gamma_e_U_.*exp(j*kzn*(z-z_)));
        I_e         =   ((Vep-Vep0).*exp(-j*kzn*(z-z_))-Vep.*Gamma_e_U_.*exp(j*kzn*(z-z_)))./Zen;
        V_h         =   ((Vhp-Vhp0).*exp(-j*kzn*(z-z_))+Vhp.*Gamma_h_U_.*exp(j*kzn*(z-z_)));
        I_h         =   ((Vhp-Vhp0).*exp(-j*kzn*(z-z_))-Vhp.*Gamma_h_U_.*exp(j*kzn*(z-z_)))./Zhn;
    end
    if z<z_
        V_e         =   ((Ven-Ven0).*exp(j*kzn*(z-z_))+Ven.*Gamma_e_D_.*exp(-j*kzn*(z-z_)));
        I_e         =   -((Ven-Ven0).*exp(j*kzn*(z-z_))-Ven.*Gamma_e_D_.*exp(-j*kzn*(z-z_)))./Zen;
        V_h         =   ((Vhn-Vhn0).*exp(j*kzn*(z-z_))+Vhn.*Gamma_h_D_.*exp(-j*kzn*(z-z_)));
        I_h         =   -((Vhn-Vhn0).*exp(j*kzn*(z-z_))-Vhn.*Gamma_h_D_.*exp(-j*kzn*(z-z_)))./Zhn;
    end
end
%% Case n<m
if n<m
    [~,Theta_n,~,~]=kpParameters(kp,m-1,Sheet);
    [Gamma_e_U,Gamma_h_U]=Refl(kp,'U',m-1,Sheet);
    zm          = 	Data(m,1);
    Vp_e      	=   Vep.*(exp(-j*kzn*(zm-z_))+Gamma_e_U_.*exp(j*kzn*(zm-z_)))./(1+Gamma_e_U.*exp(-j*2*Theta_n));
    Vp_h      	=   Vhp.*(exp(-j*kzn*(zm-z_))+Gamma_h_U_.*exp(j*kzn*(zm-z_)))./(1+Gamma_h_U.*exp(-j*2*Theta_n));
    [kzn,Theta_n,Zen,Zhn]=kpParameters(kp,n,Sheet);
    zm          = 	Data(n,1);
    [Gamma_e_U,Gamma_h_U]=Refl(kp,'U',n,Sheet);
    [Vn_e,Vn_h]=tau(kp,Vp_e,Vp_h,'U',m,n,Sheet);
    V_e         =   Vn_e.*exp(-j*Theta_n).*(exp(-j*kzn*(z-zm))+Gamma_e_U.*exp(j*kzn*(z-zm)));  
    I_e         =   Vn_e.*exp(-j*Theta_n).*(exp(-j*kzn*(z-zm))-Gamma_e_U.*exp(j*kzn*(z-zm)))./Zen;
    V_h         =   Vn_h.*exp(-j*Theta_n).*(exp(-j*kzn*(z-zm))+Gamma_h_U.*exp(j*kzn*(z-zm)));  
    I_h         =   Vn_h.*exp(-j*Theta_n).*(exp(-j*kzn*(z-zm))-Gamma_h_U.*exp(j*kzn*(z-zm)))./Zhn;
end
%% Case n>m
if n>m
    [~,Theta_n,~,~]=kpParameters(kp,m+1,Sheet);
    [Gamma_e_D,Gamma_h_D]=Refl(kp,'D',m+1,Sheet);
    zn          = 	Data(m+1,1);
    Vn_e      	=   Ven.*(exp(j*kzn*(zn-z_))+Gamma_e_D_.*exp(-j*kzn*(zn-z_)))./(1+Gamma_e_D.*exp(-j*2*Theta_n));
    Vn_h      	=   Vhn.*(exp(j*kzn*(zn-z_))+Gamma_h_D_.*exp(-j*kzn*(zn-z_)))./(1+Gamma_h_D.*exp(-j*2*Theta_n));
    [kzn,Theta_n,Zen,Zhn]=kpParameters(kp,n,Sheet);
    zn          = 	Data(n+1,1);
    [Gamma_e_D,Gamma_h_D]=Refl(kp,'D',n,Sheet);
    [Vn_e,Vn_h]=tau(kp,Vn_e,Vn_h,'D',m,n,Sheet);
    V_e         =   Vn_e.*exp(-j*Theta_n).*(exp(j*kzn*(z-zn))+Gamma_e_D.*exp(-j*kzn*(z-zn)));  
    I_e         =   -Vn_e.*exp(-j*Theta_n).*(exp(j*kzn*(z-zn))-Gamma_e_D.*exp(-j*kzn*(z-zn)))./Zen;
    V_h         =   Vn_h.*exp(-j*Theta_n).*(exp(j*kzn*(z-zn))+Gamma_h_D.*exp(-j*kzn*(z-zn)));  
    I_h         =   -Vn_h.*exp(-j*Theta_n).*(exp(j*kzn*(z-zn))-Gamma_h_D.*exp(-j*kzn*(z-zn)))./Zhn;
end
end
