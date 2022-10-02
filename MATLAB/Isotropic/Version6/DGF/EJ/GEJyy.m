function[I]=GEJyy(rho,phi,z,z_)
I1      =   -0.5*GEJ1(rho,z,z_,0);
I2      =   -0.5*cos(2*phi)*GEJ2(rho,z,z_,2);
I       =   I1+I2;
%% Add Free-Space Solution
m       =   SelectLayer(z_);
n       =   SelectLayer(z);
if m==n
    x       =   rho*cos(phi);
    y       =   rho*sin(phi);
    [~,~,~,~,~,~,Data]=Configs();
    k       =   Data(m+1,4);
    mu      =   Data(m+1,2);
    [~,~,~,~,Gyy,~,~,~,~]=GEJ0(x,y,z,z_,k,mu);
    I       =   Gyy+I;  
end
end
%%