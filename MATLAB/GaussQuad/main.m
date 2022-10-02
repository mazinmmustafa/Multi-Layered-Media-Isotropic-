close all; clear; clc
%%
n       =   16;
%% Gauss-Legendre
[x,w]	=	GaussLegendre(n);
%%
fprintf("x = \n");
for i=1:n
    if (x(i)<0)
        fprintf("-%1.14e, ", -x(i));
    else
        fprintf("+%1.14e, ", +x(i));
    end
    if (mod(i,3)==0)
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
    if (mod(i,3)==0)
        fprintf("\n");
%         fprintf("...\n");
    end
end
fprintf("\n\n");
%%