function[r1,r2,r3,r4]=CIM(pol,positions,N_max,N_boxes_max)
switch nargin
    case 2
        N_max       = 10;
        N_boxes_max = 200;
	case 3
        N_boxes_max = 200;
end
j                   =   sqrt(-1);
fprintf('Meshing and generating boxes...\n');
[Boxes]             =   ContourParameters(pol,positions,N_max,N_boxes_max);
[N_boxes,~]         =   size(Boxes);
clc;
%%
r                   =   [];
fprintf('Solving for roots...\n');
for i=1:N_boxes
    %%
    x1              =   Boxes(i,2);
    y1              =   Boxes(i,3);
    x2              =   Boxes(i,4);
    y2              =   Boxes(i,5);
    positions       =   [x1,y1,x2,y2];
	%%
    fprintf('Step:\t%0.0f/100\n',i*100/N_boxes);
	S               =   DelvesLyness(pol,positions);
    S               =   transpose(S);
    r               =   horzcat(r,S);
end
%%
clc;
r1              =   zeros(length(r));
r2              =   zeros(length(r));
r3              =   zeros(length(r));
r4              =   zeros(length(r));
for i=1:length(r)
    fprintf('Polishing the resutls...\n');
    r1(i)         	=   NewtonRaphson(@(kp)DFRiemann(kp,pol,1),adiff(@(kp)DFRiemann(kp,pol,1),1),r(i));
    fprintf('Polishing the resutls...\n');
	r2(i)         	=   NewtonRaphson(@(kp)DFRiemann(kp,pol,2),adiff(@(kp)DFRiemann(kp,pol,2),1),r(i));
    fprintf('Polishing the resutls...\n');
    r3(i)         	=   NewtonRaphson(@(kp)DFRiemann(kp,pol,3),adiff(@(kp)DFRiemann(kp,pol,3),1),r(i));
    fprintf('Polishing the resutls...\n');
    r4(i)         	=   NewtonRaphson(@(kp)DFRiemann(kp,pol,4),adiff(@(kp)DFRiemann(kp,pol,4),1),r(i));
end
NN                  =   14;
for i=1:length(r1)
    N_digits        =   NN-round(log10(abs(r1(i))));
    N_digits(isinf(N_digits))=14;
    r1_Re        	=   round(real(r1(i)),N_digits);
    r1_Im          	=   round(imag(r1(i)),N_digits);
    r1(i)          	=   r1_Re+j*r1_Im;
end
for i=1:length(r2)
    N_digits        =   NN-round(log10(abs(r2(i))));
    N_digits(isinf(N_digits))=14;
    r2_Re        	=   round(real(r2(i)),N_digits);
    r2_Im          	=   round(imag(r2(i)),N_digits);
    r2(i)          	=   r2_Re+j*r2_Im;
end
for i=1:length(r3)
    N_digits        =   NN-round(log10(abs(r3(i))));
    N_digits(isinf(N_digits))=14;
    r3_Re        	=   round(real(r2(i)),N_digits);
    r3_Im          	=   round(imag(r2(i)),N_digits);
    r3(i)          	=   r3_Re+j*r3_Im;
end
for i=1:length(r4)
    N_digits        =   NN-round(log10(abs(r4(i))));
    N_digits(isinf(N_digits))=14;
    r4_Re        	=   round(real(r4(i)),N_digits);
    r4_Im          	=   round(imag(r4(i)),N_digits);
    r4(i)          	=   r4_Re+j*r4_Im;
end
r1(abs(r1)<1e-16) 	=   0;
r2(abs(r2)<1e-16) 	=   0;
r3(abs(r3)<1e-16) 	=   0;
r4(abs(r4)<1e-16) 	=   0;
r1               	=   RemoveDuplicates(r1);
r2               	=   RemoveDuplicates(r2);
r3               	=   RemoveDuplicates(r3);
r4               	=   RemoveDuplicates(r4);
r1                 	=   sort(r1,'ComparisonMethod','real');
r2                 	=   sort(r2,'ComparisonMethod','real');
r3                 	=   sort(r3,'ComparisonMethod','real');
r4                 	=   sort(r4,'ComparisonMethod','real');
r1                 	=   RemoveZeros(r1);
r2                 	=   RemoveZeros(r2);
r3                 	=   RemoveZeros(r3);
r4                 	=   RemoveZeros(r4);
end
function[Boxes]=ContourParameters(pol,positions,N_max,N_boxes_max)
dlta            =   1e-3;
N_trials        =   2;
Boxes           =   zeros(1,5);
%%
x1              =   positions(1);
y1              =   positions(2);
x2              =   positions(3);
y2              =   positions(4);
%%
Boxes(1,1)      =   1;
Boxes(1,2)      =   x1;
Boxes(1,3)      =   y1;
Boxes(1,4)      =   x2;
Boxes(1,5)      =   y2;
[N_boxes,~]     =   size(Boxes);
count           =   0;
while (N_boxes<=N_boxes_max && count<N_boxes)
    i        	=   1;
    while (i<=N_boxes)
        [N_boxes,~]             =   size(Boxes);
        if N_boxes>N_boxes_max
            error('Maximum number of boxes exceeded!');
        end
        count_trials            =   0;
        while (count_trials<=N_trials)
            x1               	=   Boxes(i,2);  
            y1                	=   Boxes(i,3);
            x2                	=   Boxes(i,4);
            y2                	=   Boxes(i,5);
            warning('off','all');
            N                 	=   mu_k(pol,0,Boxes(i,2:5));
            N(isnan(N))         =   0.5;
            warning('on','all');
            fprintf('N\t=\t%0.4f\t+j\t%0.4f\n',real(N),imag(N));
            if ~mod(round(real(N),4),1) && (round(imag(N),4)==0) 
                Boxes(i,1)    	=   round(N);
                break;
            else
                Boxes(i,2)    	=   x1-dlta;
                Boxes(i,3)    	=   y1-dlta;
                Boxes(i,4)    	=   x2+dlta;
                Boxes(i,5)    	=   y2+dlta; 
            end
            count_trials     	=   count_trials+1;
            if count_trials>N_trials
                error('Maximum number of trials exceeded!');
            end
        end
        if round(N)==0 || round(N)<=N_max
            Boxes(i,1)          =   round(N);
            count               =   count+1;
        end
        if round(N)>N_max
            Boxes_temp          =   zeros(4,5);
            %% New box 1
            Boxes_temp(1,2)     =   x1;
            Boxes_temp(1,3)     =   y1;
            Boxes_temp(1,4)     =   x1+(x2-x1)/2;
            Boxes_temp(1,5)     =   y1+(y2-y1)/2;
            %% New box 2
            Boxes_temp(2,2)     =   x1+(x2-x1)/2;
            Boxes_temp(2,3)     =   y1;
            Boxes_temp(2,4)     =   x2;
            Boxes_temp(2,5)     =   y1+(y2-y1)/2;
            %% New box 3
            Boxes_temp(3,2)     =   x1;
            Boxes_temp(3,3)     =   y1+(y2-y1)/2;
            Boxes_temp(3,4)     =   x1+(x2-x1)/2;
            Boxes_temp(3,5)     =   y2;
            %% New box 4
            Boxes_temp(4,2)     =   x1+(x2-x1)/2;
            Boxes_temp(4,3)     =   y1+(y2-y1)/2;
            Boxes_temp(4,4)     =   x2;
            Boxes_temp(4,5)     =   y2;
            %% Insert the new boxes
            if i>1 && N_boxes>1
                Boxes           =   vertcat(Boxes(1:i-1,:),Boxes_temp,Boxes(i+1:N_boxes,:));
            end
            if i==N_boxes && N_boxes>1
                Boxes           =   vertcat(Boxes(1:i-1,:),Boxes_temp);
            end
            if i==1 && N_boxes>1
                Boxes           =   vertcat(Boxes_temp,Boxes(i+1:N_boxes,:));
            end
            if i==1 && N_boxes==1
                Boxes           =   Boxes_temp;
            end
            i                   =   i-1;
        end
        i                       =   i+1;
    end
end
end
%%
function[r]=DelvesLyness(pol,positions)
N               =   mu_k(pol,0,positions);
N               =   round(N,4);
mu_k_           =   zeros(1,N+1);
mu_k_(1)    	=   N;      
for k=1:N
    mu_k_(k+1)  =   mu_k(pol,k,positions);
end
sigma           =   zeros(1,N+1);
sigma(1)        =   1;
for k=1:N
    sum         =   0;
    for j=1:k-1
        sum     =   sum+sigma(k+1-j)*mu_k_(j+1);
    end
    sigma(k+1)  =   (-1/k)*( mu_k_(k+1) + sum );
end
r               =   roots(sigma);
end
function[I]=mu_k(pol,k,positions)
%%
j               =   sqrt(-1);
x1              =   positions(1);
y1              =   positions(2);
x2              =   positions(3);
y2              =   positions(4);
%% Sheet I
f           	=   @(kp)DFRiemann(kp,pol,1);
f_            	=   adiff(@(kp)DFRiemann(kp,pol,1),1);
IS1           	=   IRiemann();
%% Sheet II
f           	=   @(kp)DFRiemann(kp,pol,2);
f_            	=   adiff(@(kp)DFRiemann(kp,pol,2),1);
IS2           	=   IRiemann();
%% Sheet III
f           	=   @(kp)DFRiemann(kp,pol,3);
f_            	=   adiff(@(kp)DFRiemann(kp,pol,3),1);
IS3           	=   IRiemann();
%% Sheet IV
f           	=   @(kp)DFRiemann(kp,pol,4);
f_            	=   adiff(@(kp)DFRiemann(kp,pol,4),1);
IS4           	=   IRiemann();
%%
I               =   IS1+IS2+IS3+IS4;
%%
function[I]=IRiemann()
tol             =   1e-4;
D1            	=   @(t)( ((t+j*y1).^k).*f_(t+j*y1)./f(t+j*y1) );
I1              =   quadgk(D1,x1,x2,'RelTol',tol);
D2            	=   @(t)( j*((x2+j*t).^k).*f_(x2+j*t)./f(x2+j*t) );
I2              =   quadgk(D2,y1,y2,'RelTol',tol);
D3            	=   @(t)( ((t+j*y2).^k).*f_(t+j*y2)./f(t+j*y2) );
I3              =   quadgk(D3,x1,x2,'RelTol',tol);
D4            	=   @(t)( j*((x1+j*t).^k).*f_(x1+j*t)./f(x1+j*t) );
I4              =   quadgk(D4,y1,y2,'RelTol',tol);
I               =   I1+I2-I3-I4;
I           	=   I/(2*pi*j);
end
end
%%
function[zn]=NewtonRaphson(f,f_,zp)
k_max           =   12;
tol             =   1e-14;
count           =   1;
while(count<k_max)
    zn          =   zp-(f(zp)/f_(zp));
    fprintf('%d:\t%0.15f\t+j\t%0.15f\n',count,real(zn),imag(zn));
    if (abs(zn-zp)<tol*abs(zp))
        break;
    end
    zp          =   zn;
    if count==k_max-1
        zn      =   0;
    end
    count       =   count+1;
end
clc;
end
%%
function[b]=RemoveDuplicates(a)
[b,m,~]             =   unique(a,'first');
[~,d]               =   sort(m);
b                   =   b(d);
end
%%
function[r]=RemoveZeros(r)
r(r==0)=[];
end
%%
