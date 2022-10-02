function[I]=GEJzy(rho,phi,z,z_)
[~,~,k0,~,~,~,Data]=Configs();
[~,~,~,eta0,~,~,~]=Constants();
j      	=   sqrt(-1);
eps  	=   Data(SelectLayer(z)+1,3);  
I1      =   (eta0/(j*k0*eps))*sin(phi)*GEJ4(rho,z,z_,1);
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
    [~,~,~,~,~,~,~,Gzy,~]=GEJ0(x,y,z,z_,k,mu);
    I       =   Gzy+I;  
end
end
%%