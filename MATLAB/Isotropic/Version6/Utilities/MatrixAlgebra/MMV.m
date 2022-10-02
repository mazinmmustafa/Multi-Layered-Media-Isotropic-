function[B]=MMV(A,V)
[r,c]    	=   size(A);
[r_,c_]    	=   size(V);
if (r~=r_ || c/2~=c_)
    error('Wrong matrix dimensions!');
end
r           =   r/2;
c           =   c/2;
%%
A11         =   A(1:r,1:c);   
A12         =   A(1:r,c+1:2*c);
A21         =   A(r+1:2*r,1:c);   
A22         =   A(r+1:2*r,c+1:2*c);
V11       	=   V(1:r,1:c);   
V21       	=   V(r+1:2*r,1:c);   
%%
B11         =   A11.*V11+A12.*V21;
B21         =   A21.*V11+A22.*V21;
%%
B           =   [ B11 ; B21 ];
end
%%