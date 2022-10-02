function[sigma]=Graphene(w,mu_c,gamma_c,T)
j           =   sqrt(-1);
%% Constants
[~,~,~,~,kB,h_,e]=Constants();
%% Modifications
gamma_c     =   gamma_c/h_;
%% Intraband
if T==0
    A               =   -j*e*abs(mu_c);
    B               =   pi*(h_^2)*(w-j*gamma_c);
    sigma_intra     =   A./B;
else
    A               =   -2*j*e*kB*T*log(2*cosh(mu_c/(2*kB*T)));
    B               =   pi*(h_^2)*(w-j*gamma_c);
    sigma_intra     =   A./B;
end
%% Interband
A           =  	1/2 + (1/pi)*atan(((h_*(w-j*gamma_c))-2*mu_c)/(2*kB*T));
B           = 	((h_*(w-j*gamma_c))+2*mu_c).^2;
C           =  	((h_*(w-j*gamma_c))-2*mu_c).^2 +(2*kB*T)^2;
D           =  	(j/(2*pi))*log(B./C);
sigma_inter =  	(e/(4*h_))*(A+D);
%% Total
sigma       =  	sigma_intra + sigma_inter;
end