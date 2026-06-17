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
SCP_time = [];
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
    SCP_time(i) = toc;
%     fprintf('ADMM运算时间为%.6f秒\n',SCP_time(i))
    SCP_Iter_Num_log(i) = SCP_Iter_Num;
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
Cost_value = Cost_Function(1000*xlog, 1e5*ulog);
Fuel_value = Fuel_Function(1e5*ulog);
fprintf('总低价函数值为%.6f\n',Cost_value);
fprintf('总能量函数值为%.6f\n',Fuel_value);
fprintf('最大SCP迭代次数为%d次\n',max(SCP_Iter_Num_log));
fprintf('平均SCP迭代次数为%.2f次\n',mean(SCP_Iter_Num_log));
fprintf('最大MPC求解时间%.6f秒\n',max(SCP_time));
fprintf('平均MPC求解时间%.6f秒\n',mean(SCP_time));

Plot_TinyMPCFigure();