function[Ex,Ey,Ez]=ComputeFields(x,y,z,J,Mesh)
%%
[N,~]      	=   size(Mesh);
Ex          =   0;
Ey          =   0;
Ez          =   0;
for P=1:N
	xm          =   Mesh(P,1);
 	ym          =   Mesh(P,2);
 	zm          =   Mesh(P,3);
	dm          =   Mesh(P,4);
 	dVm         =   dm^3;
 	[phi,rho]   =   cart2pol(x-xm,y-ym);  
	Jx          =   J(3*(P-1)+1);
	Jy          =   J(3*(P-1)+2);
  	Jz          =   J(3*(P-1)+3);
  	Ex          =   Ex+dVm*GEJxx(rho,phi,z,zm)*Jx+dVm*GEJxy(rho,phi,z,zm)*Jy+dVm*GEJxz(rho,phi,z,zm)*Jz;
 	Ey          =   Ey+dVm*GEJyx(rho,phi,z,zm)*Jx+dVm*GEJyy(rho,phi,z,zm)*Jy+dVm*GEJyz(rho,phi,z,zm)*Jz;
   	Ez          =   Ez+dVm*GEJzx(rho,phi,z,zm)*Jx+dVm*GEJzy(rho,phi,z,zm)*Jy+dVm*GEJzz(rho,phi,z,zm)*Jz;
end
end
%%