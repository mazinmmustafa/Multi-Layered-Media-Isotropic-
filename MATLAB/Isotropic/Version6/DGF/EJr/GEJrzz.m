function[I]=GEJrzz(rho,~,z,z_)
[~,~,k0,~,~,~,Data]=Configs();
[~,~,~,eta0,~,~,~]=Constants();
eps   	=   Data(SelectLayer(z)+1,3);  
eps_   	=   Data(SelectLayer(z_)+1,3); 
I1      =   -((eta0/k0)^2/(eps*eps_))*GEJ5(rho,z,z_,0);
I       =   I1;
end
%%