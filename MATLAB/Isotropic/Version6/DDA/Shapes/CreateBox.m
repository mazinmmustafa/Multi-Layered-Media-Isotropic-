function[Data]=CreateBox(xc,yc,zc,Lx,Ly,Lz,epsn,option)
%%
[~,lambda0,~,~,~,~,~]=Configs();
if option==1
    dn          =   lambda0/(2*pi*abs(sqrt(epsn)));
end
if option==2
    dn          =   lambda0/(3*pi*abs(sqrt(epsn)));
end
if option==3
    dn          =   lambda0/(4*pi*abs(sqrt(epsn)));
end
Nx          =   round(Lx/dn);
Ny          =   round(Ly/dn);
Nz          =   round(Lz/dn);
Lx_         =   Nx*dn;
Ly_         =   Ny*dn;
Lz_         =   Nz*dn;
N           =   Nx*Ny*Nz;
Data        =   zeros(N,5);
count       =   1;
for i=1:Nx  
    for ii=1:Ny
        for iii=1:Nz
            xn          =   xc-Lx_/2+dn/2+(i-1)*dn;
            yn          =   yc-Ly_/2+dn/2+(ii-1)*dn;
            zn          =   zc-Lz_/2+dn/2+(iii-1)*dn;
            Data(count,1)	=   xn;
            Data(count,2)	=   yn;
            Data(count,3)	=   zn;
            Data(count,4)	=   dn;
            Data(count,5)	=   epsn;
            count           =   count+1;
        end
    end
end
end
%%