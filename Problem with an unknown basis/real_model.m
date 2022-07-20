function dx = real_model(x,u)
q2=x(2);
dq1=x(3);
dq2=x(4);
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;
M = [p1+2*p3*cos(q2) p2+p3*cos(q2);
     p2+p3*cos(q2) p2];
Vm = [-p3*sin(q2)*dq2 -p3*sin(q2)*(dq1+dq2);
      p3*sin(q2)*dq1 0];
Fd = diag([5.3 1.1]);
Fs = [8.45*tanh(dq1);2.35*tanh(dq2)];
d2q = pinv(M)*(u-Vm*[dq1;dq2]-Fd*[dq1;dq2]-Fs);
dx = [dq1;dq2;d2q(1);d2q(2)];
end