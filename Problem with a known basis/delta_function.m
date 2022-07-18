function delta_mau = delta_function(x,Wc,theta_mau)
Q = eye(2);
R = 1;
[~,d_sigma] = basis_function(x);
g = g_function(x);
Y = Y_function(x);
delta_mau = x'*Q*x - 1/4*Wc'*d_sigma*g*pinv(R)*g'*d_sigma'*Wc + ...
    Wc'*d_sigma*(Y*theta_mau);
end