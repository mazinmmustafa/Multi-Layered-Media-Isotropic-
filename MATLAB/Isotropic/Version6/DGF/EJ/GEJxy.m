function[I]=GEJxy(rho,phi,z,z_)
I1      =   0.5*sin(2*phi)*GEJ2(rho,z,z_,2);
I       =   I1;
%% Add Free-Space Solution
m       =   SelectLayer(z_);
n       =   SelectLayer(z);
if m==n
    x       =   rho*cos(phi);
    y       =   rho*sin(phi);
    [~,~,~,~,~,~,Data]=Configs();
    k       =   Data(m+1,4);
    mu      =   Data(m+1,2);
    [~,Gxy,~,~,~,~,~,~,~]=GEJ0(x,y,z,z_,k,mu);
    I       =   Gxy+I;  
end
end
%%