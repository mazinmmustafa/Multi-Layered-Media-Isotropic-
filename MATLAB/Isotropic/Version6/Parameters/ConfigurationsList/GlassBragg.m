function[N,lambda0,k0,omega,Gamma_U,Gamma_D,Data]=GlassBragg()
%% Load Constants
[~,~,~,~,nm,~,~,~,~]=Units();
[c0,~,~,eta0,~,~,~]=Constants();
%% General Definitions
lambda0         =   633*nm;
k0              =   2*pi/lambda0;
omega           =   k0*c0;
N               =   19;
Gamma_U         =   0;
Gamma_D         =   0;
%%
epsL            =   2.1229;
epsH            =   4.5967;
%% Create Configuration Data
j               =   sqrt(-1);
Data            =   zeros(N+1,6);
n               =   0;
% Layer 0
n               =   n+1;
z               =   2500*nm;
mu              =   0;
eps             =   0;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 1
n               =   n+1;
z               =   1523*nm;
mu              =   1;
eps             =   2.3013;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 2
n               =   n+1;
z               =   1445*nm;
mu              =   1;
eps             =   epsH;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 3
n               =   n+1;
z               =   1329*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 4
n               =   n+1;
z               =   1241*nm;
mu              =   1;
eps             =   epsH;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 5
n               =   n+1;
z               =   1115*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 6
n               =   n+1;
z               =   1037*nm;
mu              =   1;
eps             =   epsH;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 7
n               =   n+1;
z               =   911*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 8
n               =   n+1;
z               =   833*nm;
mu              =   1;
eps             =   epsH;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 9
n               =   n+1;
z               =   707*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 10
n               =   n+1;
z               =   629*nm;
mu              =   1;
eps             =   epsH;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 11
n               =   n+1;
z               =   503*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 12
n               =   n+1;
z               =   425*nm;
mu              =   1;
eps             =   epsH;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 13
n               =   n+1;
z               =   299*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 14
n               =   n+1;
z               =   221*nm;
mu              =   1;
eps             =   epsH;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 15
n               =   n+1;
z               =   69*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 16
n               =   n+1;
z               =   42*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 17
n               =   n+1;
z               =   0*nm;
mu              =   1;
eps             =   -18.3511-j*0.4331;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 18
n               =   n+1;
z               =   -27*nm;
mu              =   1;
eps             =   epsL;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 19
n               =   n+1;
z               =   -500*nm;
mu              =   1;
eps             =   1;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
end
%%