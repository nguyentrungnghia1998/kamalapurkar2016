%% Problem with unknown basis
clc; clear ; close all
%% Time step 
Step = 0.001;
T_end = 10;
t= 0:Step:T_end;
%% Variable
x = cell(1,length(t));
Wc = cell(1,length(t));
GAMMA = cell(1,length(t));
u = cell(1,length(t));
dx = cell(1,length(t));
%% Parameter
theta_approximate = load('theta_approximate.mat');
theta_approximate = theta_approximate.theta_approximate;
p = 0;
p_ = 30;
l = 1;
R = eye(2);
Q = diag([10 10 1 1]);
eps = 0.5;
eta_c1 = 1;
eta_c2 = 1;
nuy = 0.0005;
%% Initial value
Wc{1} = [5;5;0;0;0;0;25;0;2;2];
x{1} = [1;1;0;0];
GAMMA{1} = 1000*eye(10);
Xk = zeros(4,p_);
Uk = zeros(2,p_);
dXk = zeros(4,p_);
%% Simulation
for i = 1:length(t)
    g = g_function(x{i});
    [~,d_sigma] = basis_function(x{i});
    Y = Y_function(x{i});
    bm = bm_function(x{i});
    u{i} = -1/2*pinv(R)*g'*d_sigma'*Wc{i};
    dx{i} = approximate_model(x{i},u{i},theta_approximate);
    if min(eig(GAMMA{i}))<1
        GAMMA{i} = GAMMA{1};
    end
    if p==0
        p = p+1;
        l = p;
        Xk(:,p) = x{i};
        Uk(:,p) = u{i};
        dXk(:,p) = dx{i};
    else
        Xk_old = Xk(:,1:p);
        Xk_new = [Xk_old x{i}];
        if (norm(x{i}-Xk(:,l))^2/norm(x{i})>=eps)||(rank(Xk_new)>rank(Xk_old))
            if p<p_
                p = p+1;
                l = p;
                Xk(:,p) = x{i};
                Uk(:,p) = u{i};
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
                    Uk(:,l) = u{i};
                    dXk(:,l) = dx{i};
                end
            end
        end
    end
    omega = d_sigma*(Y*theta_approximate + bm + g*u{i});
    ro = 1+nuy * omega'*GAMMA{i}*omega;
    delta_mau = delta_function(x{i},Wc{i},theta_approximate);
    dWc = -eta_c1*GAMMA{i}*omega/ro*delta_mau;
    for j = 1:p
        [~,d_sigma_i] = basis_function(Xk(:,j));
        Y_i = Y_function(Xk(:,j));
        g_i = g_function(Xk(:,j));
        bm_i = bm_function(Xk(:,j));
        u_i = -1/2*pinv(R)*g_i'*d_sigma_i'*Wc{i};
        omega_i = d_sigma_i*(Y_i*theta_approximate + bm_i +g_i*u_i);
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
plot(t,xm);