close all; clear; clc
%%
global x_ w_
%%
n       =   64;
%% Gauss-Legendre
[x,w]	=	GaussLegendre(n);
x_      =   x';
w_      =   w';
%%
fprintf("x = \n");
for i=1:n
    if (x(i)<0)
        fprintf("-%1.14e, ", -x(i));
    else
        fprintf("+%1.14e, ", +x(i));
    end
    if (mod(i,2)==0)
        fprintf("\n");
%         fprintf("...\n");
    end
end
fprintf("\n\n");
%%  
fprintf("w = \n");
for i=1:n
    if (w(i)<0)
        fprintf("-%1.14e, ", -w(i));
    else
        fprintf("+%1.14e, ", +w(i));
    end
    if (mod(i,2)==0)
        fprintf("\n");
%         fprintf("...\n");
    end
end
fprintf("\n\n");
%%
I       =   Quad(@(x)func1(x),-1.6,3.7);
fprintf('%1.14E +j %1.14E\n',real(I),imag(I));
I       =   Quad2(@(x,y)func2(x,y),0.1,2.2,-1.7,3.9);
fprintf('%1.14E +j %1.14E\n',real(I),imag(I));
%%
function[y]=func1(x)
j       =   sqrt(-1);
y       =   cos(-1.6*x.^2+j*3.7);
end
%%
function[y]=func2(x,y)
j       =   sqrt(-1);
y       =   cos(x.^2).*exp(-j*y);
end
%%
function[I]=Quad(func,a,b)
global x_ w_
hm    	=   (b-a)/2;
hp     	=   (b+a)/2;
f       =   func(hm*x_+hp);
I       =   hm*(f*w_');
end
%%
function[I]=Quad2(func,a2,b2,a1,b1)
global x_ w_
[x,y] 	=   meshgrid(x_,x_);
%%
hm1   	=   (b1-a1)/2;
hp1   	=   (b1+a1)/2;
hm2   	=   (b2-a2)/2;
hp2   	=   (b2+a2)/2;
f       =   func(hm1*x+hp1,hm2*y+hp2);
f(isnan(f))=0;
I       =   hm1*hm2*(w_*f*w_');
end
%%
