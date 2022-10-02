function[kz,Theta,Ze,Zh]=kpParameters(kp,n,Sheet)
%% Load Data
[N,~,~,omega,~,~,Data]=Configs();
[~,mu0,eps0,~,~,~,~]=Constants();
mu      =   Data(n+1,2);  
eps   	=   Data(n+1,3);
k   	=   Data(n+1,4);
kz      =   SqrtRiemann(k^2-kp.^2,1);
%% Sheet I
if Sheet==1
    if n==1
        kz      =   SqrtRiemann(k^2-kp.^2,1);   
    end
    if n==N
        kz      =   SqrtRiemann(k^2-kp.^2,1);    
    end
end
%% Sheet II
if Sheet==2
    if n==1
        kz      =   SqrtRiemann(k^2-kp.^2,2);   
    end
    if n==N
        kz      =   SqrtRiemann(k^2-kp.^2,1);    
    end
end
%% Sheet III
if Sheet==3
    if n==1
        kz      =   SqrtRiemann(k^2-kp.^2,1);   
    end
    if n==N
        kz      =   SqrtRiemann(k^2-kp.^2,2);    
    end
end
%% Sheet IV
if Sheet==4
    if n==1
        kz      =   SqrtRiemann(k^2-kp.^2,2);   
    end
    if n==N
        kz      =   SqrtRiemann(k^2-kp.^2,2);    
    end
end
%% Compute Output
Ze      =   kz./(omega*eps0*eps);
Zh      =   (omega*mu0*mu)./kz;
d      	=   Data(n,1)-Data(n+1,1); 
Theta  	=   kz.*d;
end
%% SQRT Riemann Sheets
function[w]=SqrtRiemann(z,Sheet)
w       =   sqrt(z);
if Sheet==1
w       =   w.*(imag(w)<=0)-w.*(imag(w)>0); % Proper Sheet
end
if Sheet==2
w       =   w.*(imag(w)>=0)-w.*(imag(w)<0); % Improper Sheet
end
end
%%
