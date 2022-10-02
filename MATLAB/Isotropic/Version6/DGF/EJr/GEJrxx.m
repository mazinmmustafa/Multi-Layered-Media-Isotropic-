function[I]=GEJrxx(rho,phi,z,z_)
I1      =   -0.5*GEJ1(rho,z,z_,0);
I2      =   0.5*cos(2*phi)*GEJ2(rho,z,z_,2);
I       =   I1+I2;
end
%%