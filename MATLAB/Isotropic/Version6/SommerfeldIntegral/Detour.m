function[factor,kp]=Detour(kp,rho,Distance)
%% Substitution
[~,lambda0,k0,~,~,~,~]=Configs();
kp         	=   k0*kp;
%% Detour
j          	=   sqrt(-1);
a          	=   k0*(Findkmin+1);
if rho>Distance 
    b          	=   min(k0,1/rho);
else
    b           =   k0;
end
if sqrt(rho^2+Distance^2)>15*lambda0
    b           =   0.0001*k0;
end
t          	=   kp;
x         	=   t;
y         	=   b*sin(pi*t/a).*(t<a)+0.*(t>=a);
kp        	=   x+j*y;
factor    	=   1+j*(pi*b/a)*cos(pi*t/a); 
end
%%



