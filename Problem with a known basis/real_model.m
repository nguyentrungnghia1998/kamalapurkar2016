function dx = real_model(x,u)
a = -1;
b = 1;
c = -0.5;
d = -0.5;
x1 = x(1);
x2 = x(2);
f = [x1 x2 0 0;
     0 0 x1 x2*(1-(cos(2*x1)+2)^2)]*[a;b;c;d];
g = [0;(cos(2*x1)+2)];
dx = f + g*u;
end