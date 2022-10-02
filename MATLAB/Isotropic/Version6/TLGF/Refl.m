function[Gamma_e,Gamma_h]=Refl(kp,Dir,n,Sheet)
j           =   sqrt(-1);
%% Load Data
[N,~,~,~,Gamma_U,Gamma_D,Data]=Configs();
%% Gamma Up
if Dir=='U'
    Gamma_e    	=   Gamma_U.*ones(size(kp));
    Gamma_h    	=   Gamma_U.*ones(size(kp));
    for i=2:1:n
        sigma_m    	=   Data(i,6);
        [~,Theta_m,Zem,Zhm]=kpParameters(kp,i-1,Sheet);
        [~,~,Zen,Zhn]=kpParameters(kp,i,Sheet);
        Omega_e    	=   (Zem.*Zen)./(Zem+Zen);
        Omega_h    	=   (Zhm.*Zhn)./(Zhm+Zhn);
        Gamma_mn_e 	=   (Zem-Zen)./(Zem+Zen);      
        Gamma_mn_h 	=   (Zhm-Zhn)./(Zhm+Zhn);
        A_e      	=   (Gamma_mn_e-Omega_e*sigma_m)+(1-Omega_e*sigma_m).*Gamma_e.*exp(-j*2*Theta_m);
        A_h      	=   (Gamma_mn_h-Omega_h*sigma_m)+(1-Omega_h*sigma_m).*Gamma_h.*exp(-j*2*Theta_m);
        B_e      	=   (1+Omega_e*sigma_m)+(Gamma_mn_e+Omega_e*sigma_m).*Gamma_e.*exp(-j*2*Theta_m);
        B_h      	=   (1+Omega_h*sigma_m)+(Gamma_mn_h+Omega_h*sigma_m).*Gamma_h.*exp(-j*2*Theta_m);
        Gamma_e     =   A_e./B_e;
        Gamma_h     =   A_h./B_h;
    end
end
%% Gamma Down
if Dir=='D'
    Gamma_e    	=   Gamma_D.*ones(size(kp));
    Gamma_h    	=   Gamma_D.*ones(size(kp));
    for i=N-1:-1:n
        sigma_n    	=   Data(i+1,6);
        [~,Theta_m,Zem,Zhm]=kpParameters(kp,i+1,Sheet);
        [~,~,Zen,Zhn]=kpParameters(kp,i,Sheet);
        Omega_e    	=   (Zem.*Zen)./(Zem+Zen);
        Omega_h    	=   (Zhm.*Zhn)./(Zhm+Zhn);
        Gamma_mn_e 	=   (Zem-Zen)./(Zem+Zen);      
        Gamma_mn_h 	=   (Zhm-Zhn)./(Zhm+Zhn);
        A_e      	=   (Gamma_mn_e-Omega_e*sigma_n)+(1-Omega_e*sigma_n).*Gamma_e.*exp(-j*2*Theta_m);
        A_h      	=   (Gamma_mn_h-Omega_h*sigma_n)+(1-Omega_h*sigma_n).*Gamma_h.*exp(-j*2*Theta_m);
        B_e      	=   (1+Omega_e*sigma_n)+(Gamma_mn_e+Omega_e*sigma_n).*Gamma_e.*exp(-j*2*Theta_m);
        B_h      	=   (1+Omega_h*sigma_n)+(Gamma_mn_h+Omega_h*sigma_n).*Gamma_h.*exp(-j*2*Theta_m);
        Gamma_e     =   A_e./B_e;
        Gamma_h     =   A_h./B_h;
    end
end
end
%%