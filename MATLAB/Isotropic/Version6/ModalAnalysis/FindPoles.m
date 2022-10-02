close all; clear; clc;
%% Definitions
Nmax            =   4;
pol             =   'e';
x1              =   0.5;
y1              =   -0.1;
x2              =   5;
y2              =   0.1;
%%
positions       =   [x1,y1,x2,y2];
tic;
[r1,r2,r3,r4]=CIM(pol,positions,Nmax);
toc; clc;
fprintf('Time Elapsed %0.2f minutes\n',toc/60);
%%
fprintf("Sheet I:\n");
for i=1:length(r1)
    fprintf("R%d\t:\t%0.15f\t+j\t%0.15f\n",i,real(r1(i)),imag(r1(i)));
end
fprintf("Sheet II:\n");
for i=1:length(r2)
    fprintf("R%d\t:\t%0.15f\t+j\t%0.15f\n",i,real(r2(i)),imag(r2(i)));
end
fprintf("Sheet III:\n");
for i=1:length(r3)
    fprintf("R%d\t:\t%0.15f\t+j\t%0.15f\n",i,real(r3(i)),imag(r3(i)));
end
fprintf("Sheet IV:\n");
for i=1:length(r4)
    fprintf("R%d\t:\t%0.15f\t+j\t%0.15f\n",i,real(r4(i)),imag(r4(i)));
end
%%

