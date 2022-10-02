function[I]=GEMxz(rho,phi,z,z_)
[~,~,k0,~,~,~,Data]=Configs();
[~,~,~,eta0,~,~,~]=Constants();
j      	=   sqrt(-1);
mu_  	=   Data(SelectLayer(z_)+1,2);  
I1      =   (1/(j*eta0*k0*mu_))*sin(phi)*GEM3(rho,z,z_,1);
I       =   I1;
%% Add Free-Space Solution
m       =   SelectLayer(z_);
n       =   SelectLayer(z);
if m==n
    x       =   rho*cos(phi);
    y       =   rho*sin(phi);
    [~,~,~,~,~,~,Data]=Configs();
    k       =   Data(m+1,4);
    [~,~,Gxz,~,~,~,~,~,~]=GEM0(x,y,z,z_,k);
    I       =   Gxz+I;  
end
end
%%