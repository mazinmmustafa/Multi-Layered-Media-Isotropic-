function[C]=MM(A,B)
[r,c]    	=   size(A);
[r_,c_]  	=   size(B);
if (r_~=r || c_~=c)
    error('Wrong matrix dimensions!');
end
r           =   r/2;
c           =   c/2;
%%
A11         =   A(1:r,1:c);   
A12         =   A(1:r,c+1:2*c);
A21         =   A(r+1:2*r,1:c);   
A22         =   A(r+1:2*r,c+1:2*c);
B11         =   B(1:r,1:c);   
B12         =   B(1:r,c+1:2*c);
B21         =   B(r+1:2*r,1:c);   
B22         =   B(r+1:2*r,c+1:2*c);
%%
C11         =   A11.*B11+A12.*B21;
C12         =   A11.*B12+A12.*B22;
C21         =   A21.*B11+A22.*B21;
C22         =   A21.*B12+A22.*B22;
%%
C           =   [ C11 C12 ; C21 C22 ];
end
%%