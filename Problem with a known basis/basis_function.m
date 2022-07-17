function [sigma,d_sigma] = basis_function(x)
x1 = x(1);
x2 = x(2);
sigma = [x1^2;x1*x2;x2^2];
d_sigma = [2*x1 0;
           x2 x1;
           0 2*x2];
end