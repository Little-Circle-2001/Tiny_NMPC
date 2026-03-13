close all;
clear; clc;

%% 初值与参数设置
% MPC 参数设置
N = 7;                      % 预测时域

% 设置初值
start0 = [7;7;-pi;0.9];
% start0 = [3;0;-pi;0.3];
target = [5.5;6.5;-7*pi/8;0.5];
state0 = start0 + (target - start0) * linspace(0, 1, N+1);
% x0 = [linspace(6,3,N+1);linspace(5,-1,N+1)];
input0 = [0;0];
input0 = input0 .* ones(1,N);
 
% SCP 求解器参数设置
SCP_MaxIter = 100;           % SCP最大迭代次数
SCP_Tol = 1e-2;             % SCP迭代容许收敛误差

% ADMM 求解器参数设置
rho = 0.3;                  % 惩罚项系数
eta = 1.6;                  % 加速系数
selfsigma = 1e-3;           % 反正是另一个系数
ADMM_MaxIter = 100;          % ADMM最大迭代次数
ADMM_Tol = 1e-2;            % ADMM迭代容许收敛误差

% PIPG 求解器参数设置
PIPG_MaxIter = 100;         % PIPG最大迭代次数
PIPG_Tol = 1e-2;            % PIPG容许收敛误差
omega = 500;                % PIPG超参数

%% PIPG_SCP求解器
[state_SCP_opt, input_SCP_opt, J_SCP, SCP_Iter_Num] = PIPG_Based_SCP(state0, input0, target,...
                                                                     SCP_MaxIter, SCP_Tol, ...
                                                                     eta, rho, selfsigma, ADMM_MaxIter, ADMM_Tol, ...
                                                                     PIPG_MaxIter, PIPG_Tol, omega);

Plot_Figure();