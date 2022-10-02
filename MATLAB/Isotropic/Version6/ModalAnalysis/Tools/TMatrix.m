function[DF]=TMatrix(kp,pol,Sheet)
[N,~,k0,~,Gamma_U,Gamma_D,~]=Configs();
kp          =   kp*k0;
T           =   Tn(kp,2,pol,Sheet);
for n=3:N
   T      	=   MM(Tn(kp,n,pol,Sheet),T);        
end
%% Geenerate dispersion function
j           =   sqrt(-1);
I           =   ones(size(kp));
O           =   zeros(size(kp));
if pol=='e'
    [~,Theta1,Ze1,~]=kpParameters(kp,1,Sheet);
    if (Gamma_U==0)
        TR1     =   [ Ze1 ; O ];   
    end
    if (Gamma_U==+1)
        TR1     =   [ Ze1.*exp(j*Theta1) ; Ze1.*exp(-j*Theta1) ];   
    end
    if (Gamma_U==-1)
        TR1     =   [ Ze1.*exp(j*Theta1) ; -Ze1.*exp(-j*Theta1) ];   
    end
    if (Gamma_D==0)
        TRN     =   [ I O ];   
    end
    if (Gamma_D==+1)
        TRN     =   [ I -I ];   
    end
    if (Gamma_D==-1)
        TRN     =   [ I I ];   
    end
end
if pol=='h'
    [~,Theta1,~,~]=kpParameters(kp,1,Sheet);
    [~,~,~,ZhN]=kpParameters(kp,N,Sheet);
    if (Gamma_U==0)
        TR1     =   [ I ; O ];   
    end
    if (Gamma_U==+1)
        TR1     =   [ exp(j*Theta1) ; exp(-j*Theta1) ];   
    end
    if (Gamma_U==-1)
        TR1     =   [ exp(j*Theta1) ; -exp(-j*Theta1) ];   
    end
    if (Gamma_D==0)
        TRN     =   [ I./ZhN O ];   
    end
    if (Gamma_D==+1)
        TRN     =   [ I./ZhN -I./ZhN ];   
    end
    if (Gamma_D==-1)
        TRN     =   [ I./ZhN I./ZhN ];   
    end
end
DF          =   MMHV(TRN,MMV(T,TR1));
end
%% Compute Tn
function[T]=Tn(kp,n,pol,Sheet)
j           =   sqrt(-1);
I           =   ones(size(kp));
[~,Thetan,~,~]=kpParameters(kp,n,Sheet);
if pol=='e'
    [~,~,Zn,~]=kpParameters(kp,n,Sheet);
    [~,~,Zm,~]=kpParameters(kp,n-1,Sheet);
end
if pol=='h'
    [~,~,~,Zn]=kpParameters(kp,n,Sheet);
    [~,~,~,Zm]=kpParameters(kp,n-1,Sheet); 
end
[~,~,~,~,~,~,Data]=Configs();
sigma_m     =   Data(n,6);
Qn          =   Zn./Zm;
T11         =   0.5*(I+Qn+Zn*sigma_m).*exp(j*Thetan);  
T12         =   0.5*(I-Qn+Zn*sigma_m).*exp(j*Thetan);  
T21         =   0.5*(I-Qn-Zn*sigma_m).*exp(-j*Thetan);  
T22         =   0.5*(I+Qn-Zn*sigma_m).*exp(-j*Thetan);  
T           =   [ T11 T12 ; T21 T22 ];
end
%%