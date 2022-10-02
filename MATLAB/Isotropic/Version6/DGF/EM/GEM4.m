function[I]=GEM4(rho,z,z_,v)
I = HankelTransform(@(kp)Integrand(kp,rho,z,z_,v),rho,abs(z-z_));
end
%%
function[G]=Integrand(kp,rho,z,z_,v)
[factor,kp]=Detour(kp,rho,abs(z-z_));
[~,I_e,~,~]=TLGFr(kp,'v',z,z_,1);
G           =   factor.*kp.*I_e.*besselj(v,kp*rho).*kp;   
end
%%