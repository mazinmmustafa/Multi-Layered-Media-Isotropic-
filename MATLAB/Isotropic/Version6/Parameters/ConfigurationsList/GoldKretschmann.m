function[N,lambda0,k0,omega,Gamma_U,Gamma_D,Data]=GoldKretschmann()
%% Load Constants
[~,~,~,~,nm,~,~,~,~]=Units();
[c0,~,~,eta0,~,~,~]=Constants();
%% General Definitions
j               =   sqrt(-1);
lambda0         =   633*nm;
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
z               =   2000*nm;
mu              =   0;
eps             =   0;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 1
n               =   n+1;
z               =   0*nm;
mu              =   1;
eps             =   2.3013;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 2
n               =   n+1;
z               =   -50*nm;
mu              =   1;
eps             =   -11.753-j*1.2596;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 3
n               =   n+1;
z               =   -2000*nm;
mu              =   1;
eps             =   1;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
end
%%