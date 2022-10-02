function[N,lambda0,k0,omega,Gamma_U,Gamma_D,Data]=DavidTwoLayers()
%% Load Constants
[c0,~,~,eta0,~,~,~]=Constants();
%% General Definitions
j               =   sqrt(-1);
lambda0         =   1;
k0              =   2*pi/lambda0;
omega           =   k0*c0;
N               =   2;
Gamma_U         =   0;
Gamma_D         =   0;
%% Create Configuration Data
Data            =   zeros(N+1,6);
n               =   0;
% Layer 0
n               =   n+1;
z               =   2*lambda0;
mu              =   0;
eps             =   0;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 1
n               =   n+1;
z               =   0*lambda0;
mu              =   1;
eps             =   1;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 2
n               =   n+1;
z               =   -2*lambda0;
mu              =   1;
eps             =   3-j*0.1;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
end
%%