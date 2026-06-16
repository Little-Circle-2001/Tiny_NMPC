close all;
clear; clc;

load('userReference.mat','A');dt = 0.2;N = 7;
tn = size(A,2);
start = [7;7;7;-0.0707*ones(3,1);0;-0.0948;-pi/4*3;zeros(3,1)]; % 初始状态
ref_fullstates;

% SCP 求解器参数设置
SCP_MaxIter = 20;           % SCP最大迭代次数
SCP_Tol = 5e-2;             % SCP迭代容许收敛误差

% ADMM 求解器参数设置
rho = 0.1;                  % 惩罚项系数
eta = 1.6;                  % 加速系数
selfsigma = 1e-3;           % 反正是另一个系数
ADMM_MaxIter = 100;          % ADMM最大迭代次数
rho_update_inteval = 25;
adaptive_rho_tolerance = 5;
% ADMM_Tol = 5e-2;            % ADMM迭代容许收敛误差
ADMM_abs_eps = 5e-2;
ADMM_rel_eps = 5e-2;

% PIPG 求解器参数设置
PIPG_MaxIter = 1000;         % PIPG最大迭代次数
PIPG_Tol = 1e-4;            % PIPG容许收敛误差
omega = 10;                % PIPG超参数

state0 = xr(:,1:N+1); % 初值估计
input0 = ur(:,1:N);
xlog = start;
ulog = [];
% tic;
for i = 1 : tn-N
    %% PIPG_SCP求解器
    target = xr(:,N+i);
    xref = xr(:,i:N+i);
    [state_SCP_opt, input_SCP_opt, J_SCP, SCP_Iter_Num] = PIPG_Based_SCP(state0, input0, target,...
                                                                         SCP_MaxIter, SCP_Tol, ...
                                                                         eta, rho, selfsigma, ADMM_MaxIter, ADMM_abs_eps, ADMM_rel_eps, ...
                                                                         rho_update_inteval, adaptive_rho_tolerance,...
                                                                         PIPG_MaxIter, PIPG_Tol, omega);
    u = input_SCP_opt(:,1,end);
    x = RK_dyn(state0(:,1),u,dt);
    SCP_time(i) = toc;
    fprintf('ADMM运算时间为%.6f秒\n',SCP_time(i))
    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
    if i < tn-N
        state0 = [x, state_SCP_opt(:,3:end,end), xr(:,N+i+1)];
        input0 = [input_SCP_opt(:,2:end,end), ur(:,N+i)];
    end
end
for i = tn-N+1 : tn-1
    u = input_SCP_opt(:,i-tn+N+1);
    x = RK_dyn(x,u,dt);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
end