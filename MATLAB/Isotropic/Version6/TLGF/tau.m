function[Vn_e,Vn_h]=tau(kp,Vm_e,Vm_h,Dir,m,n,Sheet)
j           =   sqrt(-1);
%%
Vn_e       	=   Vm_e;
Vn_h       	=   Vm_h;
%% tau Up
if Dir=='U'
    for i=m-2:-1:n
        [~,Theta_n,~,~]=kpParameters(kp,i,Sheet);
        [~,Theta_m,~,~]=kpParameters(kp,i+1,Sheet);
        [Gamma_n_e,Gamma_n_h]=Refl(kp,'U',i,Sheet);
        [Gamma_m_e,Gamma_m_h]=Refl(kp,'U',i+1,Sheet);
        tau_e       =   (1+Gamma_m_e).*exp(-j*Theta_m)./(1+Gamma_n_e.*exp(-j*2*Theta_n));
        tau_h       =   (1+Gamma_m_h).*exp(-j*Theta_m)./(1+Gamma_n_h.*exp(-j*2*Theta_n));
        Vn_e    	=   tau_e.*Vn_e;
        Vn_h    	=   tau_h.*Vn_h;
    end
end
%% tau Down
if Dir=='D'
    for i=m+2:1:n
        [~,Theta_n,~,~]=kpParameters(kp,i,Sheet);
        [~,Theta_m,~,~]=kpParameters(kp,i-1,Sheet);
        [Gamma_n_e,Gamma_n_h]=Refl(kp,'D',i,Sheet);
        [Gamma_m_e,Gamma_m_h]=Refl(kp,'D',i-1,Sheet);
        tau_e       =   (1+Gamma_m_e).*exp(-j*Theta_m)./(1+Gamma_n_e.*exp(-j*2*Theta_n));
        tau_h       =   (1+Gamma_m_h).*exp(-j*Theta_m)./(1+Gamma_n_h.*exp(-j*2*Theta_n));
        Vn_e    	=   tau_e.*Vn_e;
        Vn_h    	=   tau_h.*Vn_h;
    end
end
end
%%