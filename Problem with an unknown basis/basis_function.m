function [sigma,d_sigma]=basis_function(x)
x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
sigma = [x1*x3;x2*x4;x3*x2;x4*x1;x1*x2;x4*x3;x1^2;x2^2;x3^2;x4^2];
d_sigma = [x3 0 x1 0;
           0 x4 0 x2;
           0 x3 x2 0;
           x4 0 0 x1;
           x2 x1 0 0;
           0 0 x4 x3;
           2*x1 0 0 0;
           0 2*x2 0 0;
           0 0 2*x3 0;
           0 0 0 2*x4];
end