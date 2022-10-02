function[I]=GEJyz(rho,phi,z,z_)
[~,~,k0,~,~,~,Data]=Configs();
[~,~,~,eta0,~,~,~]=Constants();
j      	=   sqrt(-1);
eps_  	=   Data(SelectLayer(z_)+1,3);  
I1      =   (eta0/(j*k0*eps_))*sin(phi)*GEJ3(rho,z,z_,1);
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
    [~,~,~,~,~,Gyz,~,~,~]=GEJ0(x,y,z,z_,k,mu);
    I       =   Gyz+I;  
end
end
%%