function dx = approximate_model(x,u,theta)
x1 = x(1);
x2 = x(2);
f = [x1 x2 0 0;
     0 0 x1 x2*(1-(cos(2*x1)+2)^2)]*theta;
g = [0;(cos(2*x1)+2)];
dx = f + g*u;
end