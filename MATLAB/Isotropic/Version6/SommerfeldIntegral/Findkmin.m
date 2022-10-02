function[k_min]=Findkmin()
[~,~,k0,~,~,~,Data]=Configs();
k         	=   Data(:,4);
k_min      	=   max(abs(real(k)))/k0;       
end
%%