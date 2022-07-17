function Y = Y_function(x)
x1 = x(1);
x2 = x(2);

Y = [x1 x2 0 0;
     0 0 x1 x2*(1-(cos(2*x1)+2)^2)];
end