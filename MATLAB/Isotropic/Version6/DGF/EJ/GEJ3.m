function[I]=GEJ3(rho,z,z_,v)
I = HankelTransform(@(kp)Integrand(kp,rho,z,z_,v),rho,abs(z-z_));
end
%%
function[G]=Integrand(kp,rho,z,z_,v)
[factor,kp]=Detour(kp,rho,abs(z-z_));
[V_e,~,~,~]=TLGFr(kp,'v',z,z_,1);
G           =   factor.*kp.*V_e.*besselj(v,kp*rho).*kp;   
end
%%