function[I]=GEJrzy(rho,phi,z,z_)
[~,~,k0,~,~,~,Data]=Configs();
[~,~,~,eta0,~,~,~]=Constants();
j      	=   sqrt(-1);
eps  	=   Data(SelectLayer(z)+1,3);  
I1      =   (eta0/(j*k0*eps))*sin(phi)*GEJ4(rho,z,z_,1);
I       =   I1;
end
%%