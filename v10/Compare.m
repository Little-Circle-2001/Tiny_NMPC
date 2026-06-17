close all;
clear; clc;
global tn; %#ok<*GVMIS> 

%% 初值与参数设置
% MPC 参数设置
N = 10;                      % 预测时域
T = 600;
dt = 6;

mu = 3.986e5;
a = 6871.393;
n = sqrt(mu/a^3);

% 参考轨迹 设置初值

load('ReferenceState.mat','A');
load('ReferenceInput.mat','B');
start = [-0.5;1;0.5;0;0;0]; % -500;1000;500;0;0;0
target = [0;0;0;0;0;0];

%% Obstacle
Obs_Pos = [-0.2 -0.4; ...
            0.4  0.7; ...
            0.15  0.35];
Obs_r = [0.1; 0.1];

tn = size(A,2);

xr = A;
ur = B;
 
% SCP 求解器参数设置
SCP_MaxIter = 20;           % SCP最大迭代次数
SCP_Tol = 6e-7;             % SCP迭代容许收敛误差

% ADMM 求解器参数设置
rho = 0.8;                  % 惩罚项系数
eta = 1.6;                  % 加速系数
selfsigma = 1e-3;           % 反正是另一个系数
ADMM_MaxIter = 500;          % ADMM最大迭代次数
rho_update_inteval = 50;
adaptive_rho_tolerance = 5;
% ADMM_Tol = 5e-2;            % ADMM迭代容许收敛误差
ADMM_abs_eps = 8e-5;
ADMM_rel_eps = 8e-5;
% ADMM_abs_eps = 1e-6;
% ADMM_rel_eps = 1e-6;

% PIPG 求解器参数设置
PIPG_MaxIter = 1000;        % PIPG最大迭代次数
PIPG_Tol = 5e-7;            % PIPG容许收敛误差
omega = 0.3;                % PIPG超参数


state0 = xr(:,1:N+1);
input0 = ur(:,1:N);
xlog = start;
ulog = [];
TinyNMPC_time = [];
Mosek_Time = [];
OSQP_Time = [];
% tic;
for i = 1 : tn-N
    %% PIPG_SCP求解器    
    target = xr(:,N+i);
    [state_SCP_opt, input_SCP_opt, J_SCP, SCP_Iter_Num] = PIPG_Based_SCP(state0, input0, target,...
                                                                         Obs_Pos, Obs_r, ...
                                                                         SCP_MaxIter, SCP_Tol, ...
                                                                         eta, rho, selfsigma, ADMM_MaxIter, ADMM_abs_eps, ADMM_rel_eps, ...
                                                                         rho_update_inteval, adaptive_rho_tolerance,...
                                                                         PIPG_MaxIter, PIPG_Tol, omega);
    TinyNMPC_time(i) = toc;
%     fprintf('ADMM运算时间为%.6f秒\n',TinyNMPC_time(i))
    u = input_SCP_opt(:,1,end);
    u(1) = max(min(8e-5,u(1)),-8e-5);
    u(2) = max(min(8e-5,u(2)),-8e-5);
    u(3) = max(min(8e-5,u(3)),-8e-5);
    x = RK_dyn(state0(:,1),u,dt,n);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
    if i < tn-N
        state0 = [x, state_SCP_opt(:,3:end,end), xr(:,N+i+1)];
        input0 = [input_SCP_opt(:,2:end,end), ur(:,N+i)];
    end
end

for i = tn-N+1 : tn-1
    u = input_SCP_opt(:,i-tn+N+1,end);
    x = RK_dyn(x,u,dt,n);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
end
% Sim_time = toc;
% fprintf('总运算时间为%.6f秒\n',Sim_time)

%% Mosek求解器

state0 = xr(:,1:N+1);
input0 = ur(:,1:N);
xlog = start;
ulog = [];

for j = 1 : tn - N
    tic;
    target = xr(:,N+j);
    SCP_IterNum = 0;
    delta_J = 10;
    J0 = Cost_Function(state0,input0);
    J_SCP = J0;
    input_SCP_opt = [];
    while SCP_IterNum < SCP_MaxIter && (abs(delta_J) > SCP_Tol)
        yalmip('clear')
        [Q, q, S, R, r, A, B, gamma, D, E, c] = Make_Model(state0, input0, target, Obs_Pos, Obs_r);
        delta_x = sdpvar(6,N+1);
        delta_u = sdpvar(3,N);
        
        % Cost function
        J = 0;
        for i = 1 : N
            J = J + 1/2 * delta_x(:,i)'*Q(:,:,i)*delta_x(:,i) + q(:,i)'*delta_x(:,i);
            J = J + 1/2 * delta_u(:,i)'*R(:,:,i)*delta_u(:,i) + r(:,i)'*delta_u(:,i);
        end
        J = J + 1/2 * delta_x(:,N+1)'*Q(:,:,N+1)*delta_x(:,N+1) + q(:,N+1)'*delta_x(:,N+1);
        
        % Constraints
        Constraints = [];
        Constraints = [Constraints; delta_x(:,1)==zeros(size(state0,1),1)];
        for i = 1 : N
            Constraints = [Constraints; delta_x(:,i+1) == A(:,:,i)*delta_x(:,i) + B(:,:,i)*delta_u(:,i) + gamma(:,i+1)];
            Constraints = [Constraints; D{i}*delta_x(:,i) + E{i}*delta_u(:,i) + c{i} <=0];
        end
        Constraints = [Constraints; D{N+1}*delta_x(:,N+1) + c{N+1} <= 0];
        % Solve
    
        ops = sdpsettings('verbose',0,'solver','Mosek');
        sol = optimize(Constraints,J,ops);
        Yalmip_time = toc;

        if sol.problem ~= 1        
            delta_x_opt = value(delta_x);
            delta_u_opt = value(delta_u);
            state1 = state0 + delta_x_opt;
            input1 = input0 + delta_u_opt;
            J1 = Cost_Function(state1,input1);
            state_SCP_opt = cat(3,state_SCP_opt,state1);
            input_SCP_opt = cat(3,input_SCP_opt,input1);
            J_SCP = cat(2, J_SCP, J1);
            delta_J = J0 - J1;
            state0 = state1;
            input0 = input1;
            J0 = J1;
            SCP_IterNum = SCP_IterNum + 1;
        else
            state_SCP_opt = [];
            input_SCP_opt = [];
            J1= [];
            break;
        end

        
        
    end

    u = input_SCP_opt(:,1,end);
    x = RK_dyn(state0(:,1),u,dt,n);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
    if j < tn-N
        state0 = [x, state_SCP_opt(:,3:end,end), xr(:,N+j+1)];
        input0 = [input_SCP_opt(:,2:end,end), ur(:,N+j)];
    end
    Mosek_Time(j) = toc;
    % fprintf('Yalmip运算时间为%.6f秒\n',MosekTime(j));

end
for j = tn-N+1 : tn-1
    u = input_SCP_opt(:,j-tn+N+1,end);
    x = RK_dyn(x,u,dt,n);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
end

%% OSQP求解器

state0 = xr(:,1:N+1);
input0 = ur(:,1:N);
xlog = start;
ulog = [];

for j = tn-N+1 : tn-1
    u = input_SCP_opt(:,j-tn+N+1,end);
    x = RK_dyn(x,u,dt,n);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
end

for j = 1 : tn - N
    tic;
    target = xr(:,N+j);
    SCP_IterNum = 0;
    delta_J = 10;
    J0 = Cost_Function(state0,input0);
    J_SCP = J0;
    input_SCP_opt = [];
    while SCP_IterNum < SCP_MaxIter && (abs(delta_J) > SCP_Tol)
        yalmip('clear')
        [Q, q, S, R, r, A, B, gamma, D, E, c] = Make_Model(state0, input0, target, Obs_Pos, Obs_r);
        delta_x = sdpvar(6,N+1);
        delta_u = sdpvar(3,N);
        
        % Cost function
        J = 0;
        for i = 1 : N
            J = J + 1/2 * delta_x(:,i)'*Q(:,:,i)*delta_x(:,i) + q(:,i)'*delta_x(:,i);
            J = J + 1/2 * delta_u(:,i)'*R(:,:,i)*delta_u(:,i) + r(:,i)'*delta_u(:,i);
        end
        J = J + 1/2 * delta_x(:,N+1)'*Q(:,:,N+1)*delta_x(:,N+1) + q(:,N+1)'*delta_x(:,N+1);
        
        % Constraints
        Constraints = [];
        Constraints = [Constraints; delta_x(:,1)==zeros(size(state0,1),1)];
        for i = 1 : N
            Constraints = [Constraints; delta_x(:,i+1) == A(:,:,i)*delta_x(:,i) + B(:,:,i)*delta_u(:,i) + gamma(:,i+1)];
            Constraints = [Constraints; D{i}*delta_x(:,i) + E{i}*delta_u(:,i) + c{i} <=0];
        end
        Constraints = [Constraints; D{N+1}*delta_x(:,N+1) + c{N+1} <= 0];
        % Solve
    
        ops = sdpsettings('verbose',0,'solver','Mosek');
        sol = optimize(Constraints,J,ops);
        Yalmip_time = toc;

        if sol.problem ~= 1        
            delta_x_opt = value(delta_x);
            delta_u_opt = value(delta_u);
            state1 = state0 + delta_x_opt;
            input1 = input0 + delta_u_opt;
            J1 = Cost_Function(state1,input1);
            state_SCP_opt = cat(3,state_SCP_opt,state1);
            input_SCP_opt = cat(3,input_SCP_opt,input1);
            J_SCP = cat(2, J_SCP, J1);
            delta_J = J0 - J1;
            state0 = state1;
            input0 = input1;
            J0 = J1;
            SCP_IterNum = SCP_IterNum + 1;
        else
            state_SCP_opt = [];
            input_SCP_opt = [];
            J1= [];
            break;
        end

        
        
    end

    u = input_SCP_opt(:,1,end);
    x = RK_dyn(state0(:,1),u,dt,n);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
    if j < tn-N
        state0 = [x, state_SCP_opt(:,3:end,end), xr(:,N+j+1)];
        input0 = [input_SCP_opt(:,2:end,end), ur(:,N+j)];
    end
    OSQP_Time(j) = toc;
    % fprintf('Yalmip运算时间为%.6f秒\n',OSQPTime(j));

end

for j = tn-N+1 : tn-1
    u = input_SCP_opt(:,j-tn+N+1,end);
    x = RK_dyn(x,u,dt,n);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
end
rng default;
figure('Renderer', 'painters', 'Position', [700 500 550 425]) %#ok<FGREN> 
% group = categorical({'Proposed'; 'Mosek'; 'OSQP'}); % 条件标签
hold on

boxchart(ones(size(TinyNMPC_time)), TinyNMPC_time, ...
    'BoxFaceColor',[0 0.5 0.5])

boxchart(2*ones(size(Mosek_Time)), Mosek_Time, ...
    'BoxFaceColor',[0.5 0 0.5])

boxchart(3*ones(size(OSQP_Time)), OSQP_Time, ...
    'BoxFaceColor',[0.5 0 0])
xticks([1 2 3]);
xticklabels({'Proposed','Mosek','OSQP'})
set(gca,'FontSize',12)
xlabel('$Solvers$','interpreter','latex','FontSize',12)
ylabel("$Solving~Time (s)$",'interpreter','latex','FontSize',12)
grid on;
save2pdf('figures/ComparedTime.pdf')