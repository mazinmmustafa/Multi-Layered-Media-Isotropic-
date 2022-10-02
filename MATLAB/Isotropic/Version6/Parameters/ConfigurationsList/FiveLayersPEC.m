function[N,lambda0,k0,omega,Gamma_U,Gamma_D,Data]=FiveLayersPEC()
%% Load Constants
[~,~,mm,~,~,~,~,GHz,~]=Units();
[c0,~,~,eta0,~,~,~]=Constants();
%% General Definitions
f               =   30*GHz;
lambda0         =   c0/f;
k0              =   2*pi/lambda0;
omega           =   k0*c0;
N               =   5;
Gamma_U         =   0;
Gamma_D         =   -1;
%%
d1              =   1.5*mm;
d2              =   1*mm;
%% Create Configuration Data
Data            =   zeros(N+1,6);
n               =   0;
% Layer 0
n               =   n+1;
z               =   0*mm;
mu              =   0;
eps             =   0;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 1
n               =   n+1;
z               =   z;
mu              =   1;
eps             =   1;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 2
n               =   n+1;
z               =   z-d1;
mu              =   1;
eps             =   2.1;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 3
n               =   n+1;
z               =   z-d2;
mu              =   1;
eps             =   12.5;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 4
n               =   n+1;
z               =   z-d1;
mu              =   1;
eps             =   9.8;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 5
n               =   n+1;
z               =   z-d1;
mu              =   1;
eps             =   8.6;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
end
%%









