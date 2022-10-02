function[I]=GEMzx(rho,phi,z,z_)
[~,~,k0,~,~,~,Data]=Configs();
[~,~,~,eta0,~,~,~]=Constants();
j      	=   sqrt(-1);
eps  	=   Data(SelectLayer(z)+1,3);  
I1      =   (j*eta0/(k0*eps))*sin(phi)*GEM4(rho,z,z_,1);
I       =   I1;
%% Add Free-Space Solution
m       =   SelectLayer(z_);
n       =   SelectLayer(z);
if m==n
    x       =   rho*cos(phi);
    y       =   rho*sin(phi);
    [~,~,~,~,~,~,Data]=Configs();
    k       =   Data(m+1,4);
    [~,~,~,~,~,~,Gzx,~,~]=GEM0(x,y,z,z_,k);
    I       =   Gzx+I;  
end
end
%%