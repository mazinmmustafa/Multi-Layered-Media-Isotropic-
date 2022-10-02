function[I]=Romberg(func,a,b,tol,kmax)
%   function[I]=Romberg(func,a,b,tol,kmax)
%   Romberg integrates function 'func' of one variable and nonsingular    
%   from 'a' to 'b' with tolerance 'tol' and maximum order of 'kmax'. 
%   0<tol<1 & kmax>0
%  	I=Romberg(@func,a,b,tol,kmax)
%   I=Romberg(@func,a,b,tol),default kmax=15
%   I=Romberg(@func,a,b),default tol=1e-10, kmax=15
%   e.g. 
%   I=Romberg(@sin,0,pi)
%   I=2.000000000000000
switch nargin
    case 3
        tol     =	1e-10;
        kmax	=	15;
	case 4
        kmax	=	15;
end
kmax    =	abs(kmax);
R       =	zeros(1,kmax+1);
err     =	1;
Ip      =	0;
R(1)	=   ((b-a)/2)*(func(a)+func(b));
k       =   1;
while(err>tol*abs(Ip))
	R(k+1)      =   trapm(func,a,b,k+1,R(k));
    for j=k:-1:1
        p       =	4^(k-j+1);
        R(j)	=	(p*R(j+1)-R(j))/(p-1);
    end
  	err 	=	abs(R(1)-Ip);
    Ip      =	R(1);
    k       =	k+1;
    if k==kmax
        warning('Maximum trials exceeded with no convergence!');
        break;
    end
end
I       =   R(1);
function[I]=trapm(func,a,b,k,Ip)
H       =	b-a;
i       =   1:2^(k-2);
I       =   func(a+(2*i-1)*H/(2^(k-1)));
I       =   sum(I);
I       =   0.5*Ip+(H/(2^(k-1)))*I;
end
end


