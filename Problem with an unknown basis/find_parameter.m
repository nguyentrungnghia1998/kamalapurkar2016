%% Problem with a unknown basis
%% Time step 
T_end = 50;
Step = 0.001;
t = 0:Step:T_end;
%% Variable
x = cell(1,length(t));
u = cell(1,length(t));
theta = cell(1,length(t));
%% Parameter
R = eye(2);
eps = 0.005;
GAMMA_theta = diag([90 50 160 50]);
k_theta = 1.5;
%% Initial value 
x{1} = [1;1;0;0];
Wa_first = [5;5;0;0;0;0;25;0;2;2];
theta{1} = [1;1;1;1];
theta{2} = theta{1};
%% Simulation
for i = 1:length(t)
    g = g_function(x{i});
    [~,d_sigma] = basis_function(x{i});
    u{i} = -1/2*pinv(R)*g'*d_sigma'*Wa_first;
    dx = real_model(x{i},u{i});
    if i == length(t)
        break
    end
    %% Update state
    x{i+1} = x{i} + Step * dx;
end

xm = cell2mat(x);
um = cell2mat(u);
plot(t,xm);
dxm = zeros(size(xm));
dxm(1,:) = xm(3,:);
dxm(2,:) = xm(4,:);
k3 = diff(xm(3,:))/0.001;
k4 = diff(xm(4,:))/0.001;
dxm(3,:) = [k3 k3(end)];
dxm(4,:) = [k4 k4(end)];
l = 1;
p = 1;
p_ = 30;
Xk = zeros(4,p_);
dXk = zeros(4,p_);
Uk = zeros(2,p_);
Xk(:,1) = xm(:,1);
dXk(:,1) = dxm(:,1);
Uk(:,1) = um(:,1);
for i = 2:length(t)
    Xk_old = Xk(:,1:p);
    Xk_new = [Xk_old xm(:,i)];
    if (norm(xm(:,i)-Xk(:,l))^2/norm(xm(:,i))>=eps)||(rank(Xk_new)>rank(Xk_old))
        if p<p_
            p = p+1;
            l = p;
            Xk(:,p) = xm(:,i);
            Uk(:,p) = um(:,i);
            dXk(:,p) = dxm(:,i);
        else
            T = Xk_old;
            S_old = min(svd(Xk_old'));
            S = zeros(p,1);
            for j = 1:p
                Xk_old(:,j) = xm(:,i);
                S(j) = min(svd(Xk_old'));
                Xk_old = T;
            end
            [S_max,l] = max(S);
            if S_max>S_old
                Xk(:,l) = xm(:,i);
                Uk(:,l) = um(:,i);
                dXk(:,l) = dxm(:,i);
            end
        end
    end
    dtheta = [0;0;0;0];
    for j = 1:p
        Y = Y_function(Xk(:,j));
        g = g_function(Xk(:,j));
        bm = bm_function(Xk(:,j));
        dtheta = dtheta + GAMMA_theta*k_theta/p*(Y'*(dXk(:,j)-g*Uk(:,j)-Y*theta{i}-bm));
    end
    if i==length(t)
        break
    end
    %% Update theta
    theta{i+1} = theta{i} + Step * dtheta;
end

theta_m = cell2mat(theta);
plot(t,theta_m);
legend('f_{d1}','f_{d2}','f_{s1}','f_{s2}');
xlabel('Time (s)');
ylabel('$$\hat{\theta}(t)$$','Interpreter','latex');
title('Unknown Parameters in Drift Dynamics');

theta_approximate = theta{end};
save('theta_approximate.mat','theta_approximate');