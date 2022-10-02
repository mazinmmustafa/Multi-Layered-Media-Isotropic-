function[eps]=NobleMetal(lambda,element)
j           =   sqrt(-1);
h           =   6.582119e-16;
eps0        =   8.854187817e-12; 
mu0         =   1.256637061e-6;
c           =   1/sqrt(mu0*eps0);
w           =   2*pi*c./lambda;
w           =   w*h; 
s           =   j*w;
if strcmp(element,'Au')
    N       =	8;
    p       =	zeros(1,N); 
    c       =   zeros(1,N);
    e       =	1; 
    d       =   5.32796e+02;
    p(1)    =   -2.06404e-01;
    p(2)    =	-3.88325e-01;
    p(3)    =	-9.04334e-01 + j*2.29151;
    p(4)    =	conj(p(3));
    p(5)    =	-1.46294e-01 + j*2.43979;
    p(6)  	=	conj(p(5));
    p(7)    =	-4.44628e-01 + j*3.970239;
    p(8)    =	conj(p(7));
    c(1)    =	-7.48794e+02;
    c(2)    =	2.21802e+02;
    c(3)    =   5.86304 + j*6.13502e-01;
    c(4)    =   conj(c(3));
    c(5)    =   1.15178e-01 + j*1.09961e-01;
    c(6)    =   conj(c(5));
    c(7)    =   3.25546e-01 - j*3.23235e-01;
    c(8)    =   conj(c(7));
end
if strcmp(element,'Ag')
    N       =   6;
    p       =   zeros(1,N); 
    c       =   zeros(1,N);
    e       =   1; 
    d       =   2.5926e3;
    p(1)    =   -3.4089e-2;
    p(2)    =   -2.4860;
    p(3)    =   -2.5434e-1 + j*3.8737;
    p(4)    =   conj(p(3));
    p(5)    =   -8.9100e-1 + j*3.9425;
    p(6)    =   conj(p(5));
    c(1)    =   -2.5956e3;
    c(2)    =   1.1961e1;
    c(3)    =   1.0284e-1 + j*3.999e-1;
    c(4)    =   conj(c(3));
    c(5)    =   3.1782 - j*5.5464e-1;
    c(6)  	=   conj(c(5));
end
sum         =   0;
for i=1:N
    sum     =   sum+c(i)./(s-p(i));
end
eps         =   e+d./s+sum;
end
