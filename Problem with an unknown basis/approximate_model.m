function dx = approximate_model(x,u,theta_mau)
Y = Y_function(x);
bm = bm_function(x);
g = g_function(x);
dx = Y*theta_mau+bm+g*u;
end