function delta = delta_function(x,Wc,theta_mau)
Q = diag([10 10 1 1]);
R = eye(2);
Y = Y_function(x);
g = g_function(x);
bm = bm_function(x);
[~,d_sigma] = basis_function(x);
delta = x'*Q*x - 1/4*Wc'*d_sigma*g*pinv(R)*g'*d_sigma'*Wc + Wc'*d_sigma*(Y*theta_mau+bm);
end