%% Problem with a known basis 
%% Time step
Step = 0.001;
T_end = 10;
t = 0:Step:T_end;
%% Variable
x = cell(1,length(t));
Wc = cell(1,length(t));
u = cell(1,length(t));
dx = cell(1,length(t));
GAMMA = cell(1,length(t));
%% Parameter
theta_approximate = load('theta_approximate.mat');
theta_approximate = theta_approximate.theta_approximate;
p_ = 30;
p = 0;
l=1;
R = 1;
Q = eye(2);
eps = 0.01;
eta_c1 = 1;
eta_c2 = 15;
nuy = 0.005;
%% Initial value
Wc{1} = [1;1;1];
x{1} = [-1;-1];
GAMMA{1} = 100*eye(3);
Xk = zeros(2,p_);
Uk = zeros(1,p_);
dXk = zeros(2,p_);
%% Simulation
for i = 1:length(t)
    g = g_function(x{i});
    [~,d_sigma] = basis_function(x{i});
    Y = Y_function(x{i});
    u{i} = -1/2*pinv(R)*g'*d_sigma'*Wc{i};
    dx{i} = approximate_model(x{i},u{i},theta_approximate);
    if p==0
        p = p+1;
        l = p;
        Xk(:,p) = x{i};
        Uk(p) = u{i};
        dXk(:,p) = dx{i};
    else
        Xk_old = Xk(:,1:p);
        Xk_new = [Xk_old x{i}];
        if (norm(x{i}-Xk(:,l))^2/norm(x{i})>=eps)||(rank(Xk_new)>rank(Xk_old))
            if p<p_
                p = p+1;
                l = p;
                Xk(:,p) = x{i};
                Uk(p) = u{i};
                dXk(:,p) = dx{i};
            else
                T = Xk_old;
                S_old = min(svd(Xk_old'));
                S = zeros(p,1);
                for j = 1:p
                    Xk_old(:,j) = x{i};
                    S(j) = min(svd(Xk_old'));
                    Xk_old = T;
                end
                [S_max,l] = max(S);
                if S_max>S_old
                    Xk(:,l) = x{i};
                    Uk(l) = u{i};
                    dXk(:,l) = dx{i};
                end
            end
        end
    end
    omega = d_sigma*(Y*theta_approximate + g*u{i});
    ro = 1+nuy * omega'*GAMMA{i}*omega;
    delta_mau = delta_function(x{i},Wc{i},theta_approximate);
    dWc = -eta_c1*GAMMA{i}*omega/ro*delta_mau;
    for j = 1:p
        [~,d_sigma_i] = basis_function(Xk(:,j));
        Y_i = Y_function(Xk(:,j));
        g_i = g_function(Xk(:,j));
        u_i = -1/2*pinv(R)*g_i'*d_sigma_i'*Wc{i};
        omega_i = d_sigma_i*(Y_i*theta_approximate+g_i*u_i);
        ro_i = 1+nuy*omega_i'*GAMMA{i}*omega_i;
        delta_mau_i = delta_function(Xk(:,j),Wc{i},theta_approximate);
        dWc = dWc - eta_c2/p*omega_i/ro_i*delta_mau_i;
    end
    dGAMMA = -eta_c1*(GAMMA{i}*(omega*omega')*GAMMA{i})/ro^2;
    if i == length(t)
        break
    end
    %% Update state, Wc, GAMMA
    x{i+1} = x{i} + Step * real_model(x{i},u{i});
    GAMMA{i+1} = GAMMA{i} + Step * dGAMMA;
    Wc{i+1} = Wc{i} + Step*dWc;
end

xm = cell2mat(x);
um = cell2mat(u);
Wcm = cell2mat(Wc);
figure(1);
plot(t,xm);
legend('x_1','x_2');
xlabel('Time (s)');
ylabel('x(t)');
title('State Trạectory');
figure(2);
plot(t,um);
xlabel('Time (s)');
ylabel('u(t)');
title('Control Trạectory');
figure(3);
plot(t,Wcm);
legend('$$\hat{W}_{c1,1}$$','$$\hat{W}_{c1,2}$$','$$\hat{W}_{c1,3}$$','Interpreter','Latex');
xlabel('Time (s)');
ylabel('$$\hat{W}_{c1}(t)$$','Interpreter','latex');
title('Value Function Weights');