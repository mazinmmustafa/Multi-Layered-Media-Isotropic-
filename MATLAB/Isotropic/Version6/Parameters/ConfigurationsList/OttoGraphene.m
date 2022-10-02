function[N,lambda0,k0,omega,Gamma_U,Gamma_D,Data]=OttoGraphene()
%% Load Constants
[~,~,~,um,~,~,~,~,THz]=Units();
[c0,~,~,eta0,~,~,~]=Constants();
%% General Definitions
j               =   sqrt(-1);
f               =   1*THz;
lambda0         =   c0/f;
k0              =   2*pi/lambda0;
omega           =   k0*c0;
N               =   3;
Gamma_U         =   0;
Gamma_D         =   0;
%% Create Configuration Data
Data            =   zeros(N+1,6);
n               =   0;
% Layer 0
n               =   n+1;
z               =   1000*um;
mu              =   0;
eps             =   0;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 1
n               =   n+1;
z               =   0*um;
mu              =   1;
eps             =   2.003^2;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 2
n               =   n+1;
z               =   -20*um;
mu              =   1;
eps             =   1;
sigma           =   0.000369059545723-j*0.015237384931248;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 3
n               =   n+1;
z               =   -1000*um;
mu              =   1;
eps             =   1.762^2;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
end
%%









