function[N,lambda0,k0,omega,Gamma_U,Gamma_D,Data]=VisserPoles()
%% Load Constants
[~,~,~,um,nm,~,~,~,~]=Units();
[c0,~,~,eta0,~,~,~]=Constants();
%% General Definitions
j               =   sqrt(-1);
lambda0         =   1.3*um;
k0              =   2*pi/lambda0;
omega           =   k0*c0;
N               =   5;
Gamma_U         =   0;
Gamma_D         =   0;
%%
d2              =   600*nm;
d3              =   400*nm;
d4              =   600*nm;
%% Create Configuration Data
Data            =   zeros(N+1,6);
n               =   0;
% Layer 0
n               =   n+1;
z               =   0*um;
mu              =   0;
eps             =   0;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 1
n               =   n+1;
z               =   z-0*um;
mu              =   1;
eps             =   1;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 2
n               =   n+1;
z               =   z-d2;
mu              =   1;
eps             =   (3.4-j*0.002)^2;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 3
n               =   n+1;
z               =   z-d3;
mu              =   1;
eps             =   (3.6+j*0.01)^2;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 4
n               =   n+1;
z               =   z-d4;
mu              =   1;
eps             =   (3.4-j*0.002)^2;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
% Layer 5
n               =   n+1;
z               =   z-0*um;
mu              =   1;
eps             =   1;
sigma           =   0;
k               =   k0*sqrt(mu)*sqrt(eps);
eta             =   eta0*sqrt(mu)/sqrt(eps);
Data(n,:)       =   [ z mu eps k eta sigma ];
end
%%









