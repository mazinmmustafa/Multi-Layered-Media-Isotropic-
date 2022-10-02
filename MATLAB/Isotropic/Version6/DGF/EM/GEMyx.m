function[I]=GEMyx(rho,phi,z,z_)
I1      =   0.5*GEM2(rho,z,z_,0);
I2      =   0.5*cos(2*phi)*GEM1(rho,z,z_,2);
I       =   I1+I2;
%% Add Free-Space Solution
m       =   SelectLayer(z_);
n       =   SelectLayer(z);
if m==n
    x       =   rho*cos(phi);
    y       =   rho*sin(phi);
    [~,~,~,~,~,~,Data]=Configs();
    k       =   Data(m+1,4);
    [~,~,~,Gyx,~,~,~,~,~]=GEM0(x,y,z,z_,k);
    I       =   Gyx+I;  
end
end
%%