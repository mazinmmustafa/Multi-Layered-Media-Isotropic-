function[B]=MMHV(H,V)
[r,c]    	=   size(V);
[r_,c_]  	=   size(H);
if (r/2~=r_ || c~=c_/2)
    error('Wrong matrix dimensions!');
end
r           =   r/2;
c           =   c_/2;
%%
V11         =   V(1:r,1:c);   
V21         =   V(r+1:2*r,1:c);   
H11         =   H(1:r,1:c);   
H12         =   H(1:r,c+1:2*c);
%%
B           =   H11.*V11+H12.*V21;
end
%%