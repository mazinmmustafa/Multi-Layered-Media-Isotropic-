function[I]=GEM3(rho,z,z_,v)
I = HankelTransform(@(kp)Integrand(kp,rho,z,z_,v),rho,abs(z-z_));
end
%%
function[G]=Integrand(kp,rho,z,z_,v)
[factor,kp]=Detour(kp,rho,abs(z-z_));
[~,~,V_h,~]=TLGFr(kp,'i',z,z_,1);
G           =   factor.*kp.*V_h.*besselj(v,kp*rho).*kp;   
end
%%