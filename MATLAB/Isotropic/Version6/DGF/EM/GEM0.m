function[Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz]=GEM0(x,y,z,z_,k)
j               =   sqrt(-1);
%%
R               =   sqrt(x.^2+y.^2+(z-z_).^2);
g               =   exp(-j*k*R)./(4*pi*R);
g1              =   -(1+j*k*R).*g./R;
%% Exract Singularity
if R==0
    Gxx             =   0;
    Gxy             =   0;
    Gxz             =   0;
    Gyx             =   0;
    Gyy             =   0;
    Gyz             =   0;
    Gzx             =   0;
    Gzy             =   0;
    Gzz             =   0;
else
    Gxx             =   0;
    Gxy             =   ((z-z_)./R).*g1;
    Gxz             =   (-y./R).*g1;
    Gyx             =   (-(z-z_)./R).*g1;
    Gyy             =   0;
    Gyz             =   (x./R).*g1;
    Gzx             =   (y./R).*g1;
    Gzy             =   (-x./R).*g1;
    Gzz             =   0;
end
end