function[f]=DF(kp,pol)
f=TMatrix(kp,pol,1).*TMatrix(kp,pol,2).*TMatrix(kp,pol,3).*TMatrix(kp,pol,4);
end
%%