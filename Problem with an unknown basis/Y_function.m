function Y = Y_function(x)
q2 = x(2);
dq1 = x(3);
dq2 = x(4);
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;
M = [p1+2*p3*cos(q2) p2+p3*cos(q2);
     p2+p3*cos(q2) p2];
invM = pinv(M);
Y = [0 0 0 0;
     0 0 0 0;
     -invM(1,1)*dq1 -invM(1,2)*dq2 -invM(1,1)*tanh(dq1) -invM(1,2)*tanh(dq2);
     -invM(2,1)*dq1 -invM(2,2)*dq2 -invM(2,1)*tanh(dq1) -invM(2,2)*tanh(dq2)];
end