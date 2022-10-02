function[I]=GEJryz(rho,phi,z,z_)
[~,~,k0,~,~,~,Data]=Configs();
[~,~,~,eta0,~,~,~]=Constants();
j      	=   sqrt(-1);
eps_  	=   Data(SelectLayer(z_)+1,3);  
I1      =   (eta0/(j*k0*eps_))*sin(phi)*GEJ3(rho,z,z_,1);
I       =   I1;
end
%%