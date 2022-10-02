function[I]=GEMyy(rho,phi,z,z_)
I1      =   0.5*sin(2*phi)*GEM1(rho,z,z_,2);
I       =   I1;
%% Add Free-Space Solution
m       =   SelectLayer(z_);
n       =   SelectLayer(z);
if m==n
    x       =   rho*cos(phi);
    y       =   rho*sin(phi);
    [~,~,~,~,~,~,Data]=Configs();
    k       =   Data(m+1,4);
    [~,~,~,~,Gyy,~,~,~,~]=GEM0(x,y,z,z_,k);
    I       =   Gyy+I;  
end
end
%%