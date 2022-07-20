function bm = bm_function(x)
q2 = x(2);
dq1 = x(3);
dq2 = x(4);
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;
M = [p1+2*p3*cos(q2) p2+p3*cos(q2);
     p2+p3*cos(q2) p2];
Vm = [-p3*sin(q2)*dq2 -p3*sin(q2)*(dq1+dq2);
      p3*sin(q2)*dq1 0];
bm = [dq1;dq2;-pinv(M)*Vm*[dq1;dq2]];
end