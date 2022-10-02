function[I]=GEJzz(rho,phi,z,z_)
[~,~,k0,~,~,~,Data]=Configs();
[~,~,~,eta0,~,~,~]=Constants();
eps   	=   Data(SelectLayer(z)+1,3);  
eps_   	=   Data(SelectLayer(z_)+1,3); 
I1      =   -((eta0/k0)^2/(eps*eps_))*GEJ5(rho,z,z_,0);
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
    [~,~,~,~,~,~,~,~,Gzz]=GEJ0(x,y,z,z_,k,mu);
    I       =   Gzz+I;  
    if (rho==0 && z==z_)
        j       =   sqrt(-1);
        I       =   I-(eta0/(j*k0*eps));
    end
end
end
%%