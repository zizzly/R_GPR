
function y = vwDE(x,U,Q,cd,lambda)

%cd = 1.3e-3;    %drag coeff
%lambda = 4e-5;              % normalized surface feedback parameter from vdW
%Q = 1.5e-5;                 % normalized net radiation parameter from vdW


Rb=x./(U.^2); 
y=Q - lambda*x - cd*U*fm(Rb)*x; 

end

