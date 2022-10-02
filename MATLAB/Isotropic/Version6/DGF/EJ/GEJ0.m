function[Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz]=GEJ0(x,y,z,z_,k,mu)
j               =   sqrt(-1);
[~,mu0,~,~,~,~,~]=Constants();
[~,~,~,omega,~,~,~]=Configs();
%%
R               =   sqrt(x.^2+y.^2+(z-z_).^2);
g               =   exp(-j*k*R)./(4*pi*R);
g1              =   -(1+j*k*R).*g./R;
g2              =   -(1+j*k*R).*g1./R+g./(R.^2);
%%
factor          =   -j*omega*mu0*mu;
%% Exract Singularity
if R==0
    Gxx             =   factor*(-1/(3*k^2));
    Gxy             =   0;
    Gxz             =   0;
    Gyx             =   0;
    Gyy             =   factor*(-1/(3*k^2));
    Gyz             =   0;
    Gzx             =   0;
    Gzy             =   0;
    Gzz             =   factor*(-1/(3*k^2));
else
    Gxx             =   factor*(g+(1/k^2)*(g1./R+(x./R).*(x./R).*(g2-g1./R)));
    Gxy             =   factor*((1/k^2)*((x./R).*(y./R).*(g2-g1./R)));
    Gxz             =   factor*((1/k^2)*((x./R).*((z-z_)./R).*(g2-g1./R)));
    Gyx             =   factor*((1/k^2)*((y./R).*(x./R).*(g2-g1./R)));
    Gyy             =   factor*(g+(1/k^2)*(g1./R+(y./R).*(y./R).*(g2-g1./R)));
    Gyz             =   factor*((1/k^2)*((y./R).*((z-z_)./R).*(g2-g1./R)));
    Gzx             =   factor*((1/k^2)*(((z-z_)./R).*(x./R).*(g2-g1./R)));
    Gzy             =   factor*((1/k^2)*(((z-z_)./R).*(y./R).*(g2-g1./R)));
    Gzz             =   factor*(g+(1/k^2)*(g1./R+((z-z_)./R).*((z-z_)./R).*(g2-g1./R)));
end
end