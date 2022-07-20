function g = g_function(x)
q2 = x(2);
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;
M = [p1+2*p3*cos(q2) p2+p3*cos(q2);
     p2+p3*cos(q2) p2];
g = [zeros(2);pinv(M)];
end