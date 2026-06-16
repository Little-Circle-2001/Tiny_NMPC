close all;
clear; clc;
global tn; %#ok<*GVMIS> 

%% 初值与参数设置
% MPC 参数设置
N = 7;                      % 预测时域
T = 20;
dt = 0.2;

%% Obstacle
Obs_Pos = [-0.200 -0.400; ...
            0.400  0.700; ...
            0.150  0.350];
Obs_r = [0.100; 0.100];

% Obs_Pos = [-200; ...
%             400; ...
%             150];
% Obs_r = [100];

%% 参考轨迹 设置初值
load('userReference.mat','A');
start = [7;7;7;-0.0707*ones(3,1);0;-0.0948;-pi/4*3;zeros(3,1)]; % 初始状态
ref_fullstates;

% xr = start;
tn = size(A,2);
% for i = 2:tn-1
%     vx = (A(1,i)-A(1,i-1))/dt;
%     vy = (A(2,i)-A(2,i-1))/dt;
%     vz = (A(3,i)-A(3,i-1))/dt;
%     xr = [xr, [A(:,i);vx;vy;vz]];
% end
% xr = [xr, target];
% ur = [];
% for i = 1:tn-1
%     ux = (xr(4,i+1)-xr(4,i))/dt - 3*n^2*xr(1,i) - 2*n*xr(5,i);
%     uy = (xr(4,i+1)-xr(4,i))/dt + 2*n*xr(4,i);
%     uz = (xr(4,i+1)-xr(4,i))/dt + n^2*xr(3,i);
%     ur = [ur, [ux;uy;uz]];
% end

state0 = xr(:,1:N+1);
input0 = ur(:,1:N);
state_SCP_opt = state0;
input_SCP_opt = [];
xlog = start;
ulog = [];

% SCP 求解器参数设置
SCP_MaxIter = 20;          % SCP最大迭代次数
SCP_Tol = 6e-7;             % SCP迭代容许收敛误差



%% SCP求解器


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
        delta_x = sdpvar(size(xr,1),N+1);
        delta_u = sdpvar(size(ur,1),N);
        
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
    
        ops = sdpsettings('verbose',0,'solver','OSQP');
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
%     u(1) = max(min(8e-5,u(1)),-8e-5);
%     u(2) = max(min(8e-5,u(2)),-8e-5);
%     u(3) = max(min(8e-5,u(3)),-8e-5);
    x = RK_dyn(state0(:,1),u,dt);

    ulog = [ulog, u]; %#ok<*AGROW> 
    xlog = [xlog, x];
    if j < tn-N
        state0 = [x, state_SCP_opt(:,3:end,end), xr(:,N+j+1)];
        input0 = [input_SCP_opt(:,2:end,end), ur(:,N+j)];
    end
SCP_time(j) = toc;
SCP_Iter_Num_log(j) = SCP_IterNum;
% fprintf('Yalmip运算时间为%.6f秒\n',SCP_time(j));

end

for j = tn-N+1 : tn-1
    u = input_SCP_opt(:,j-tn+N+1,end);
    x = RK_dyn(x,u,dt);

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

figure('Renderer', 'painters', 'Position', [600 100 800 700]) %#ok<FGREN> 

% plot obstacle
for k = 1 : length(Obs_r)
    [obs_x, obs_y, obs_z] = sphere(50); % 50是网格的分辨率
    obs_x = 1000 * Obs_r(k) * obs_x + 1000 * Obs_Pos(1,k);
    obs_y = 1000 * Obs_r(k) * obs_y + 1000 * Obs_Pos(2,k);
    obs_z = 1000 * Obs_r(k) * obs_z + 1000 * Obs_Pos(3,k);
    
    surf(obs_x, obs_y, obs_z); % 使用'none'去除边缘线
    colormap(hot); % 使用hot颜色映射
    
    axis equal; % 保持轴比例相等，使球体看起来是圆形的
    
    hold on;
end

h1 = plot3(1000 * xlog(1,:),1000 * xlog(2,:),1000 * xlog(3,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r');
h2 = plot3(1000 * xr(1,:),1000 * xr(2,:),1000 * xr(3,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b');
grid on;
xlabel("$x_1 (m)$",'interpreter','latex','FontSize',12)
ylabel("$x_2 (m)$",'interpreter','latex','FontSize',12)
legend([h1,h2],"Actual","Reference",'interpreter','latex','FontSize',12, 'Location', 'northwest')
view(48,34)
set(gca, 'FontSize', 12);
save2pdf('figures/MosekTrajectory.pdf')

figure('Renderer', 'painters', 'Position', [700 500 1100 425]) %#ok<FGREN>
subplot(1,2,1)
plot(1000 * xlog(4,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
hold on;
plot(1000 * xlog(5,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','g')
plot(1000 * xlog(6,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','b')
% plot(xr(4,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b');
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$v (m/s)$",'interpreter','latex','FontSize',12)
legend("$v_x$","$v_y$","$v_z$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);

subplot(1,2,2)
plot(1000 * xr(4,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','r');
hold on;
plot(1000 * xr(5,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','g');
plot(1000 * xr(6,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b');
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$v_r (m/s)$",'interpreter','latex','FontSize',12)
legend("$v_{rx}$","$v_{ry}$","$v_{rz}$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);
save2pdf('figures/MosekVelocity.pdf')


figure('Renderer', 'painters', 'Position', [700 500 1100 425]) %#ok<FGREN> 
subplot(1,2,1)
hold on;
plot(1e5 * ulog(1,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
plot(1e5 * ulog(2,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','g')
plot(1e5 * ulog(3,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','b')
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$u (m/s^2)$",'interpreter','latex','FontSize',12)
legend("$u_x$","$u_y$","$u_z$",'interpreter','latex','FontSize',12);

subplot(1,2,2)
hold on;
plot(1e5 * ur(1,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','r')
plot(1e5 * ur(2,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','g')
plot(1e5 * ur(3,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b')
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$u_r (m/s^2)$",'interpreter','latex','FontSize',12)
legend("$u_{rx}$","$u_{ry}$","$u_{rz}$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);
save2pdf('figures/MosekInput.pdf')

% figure('Renderer', 'painters', 'Position', [700 500 550 425]) %#ok<FGREN> 
% hold on;
% plot(ulog(2,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
% plot(ur(2,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b')
% grid on;
% xlabel("$t (s)$",'interpreter','latex','FontSize',12)
% ylabel("$a (m \cdot s^{-2})$",'interpreter','latex','FontSize',12)
% legend("$a$",'interpreter','latex','FontSize',12);
% set(gca, 'FontSize', 12);

% figure('Renderer', 'painters', 'Position', [200 100 1200 600]) %#ok<FGREN> 
% plot(0:SCP_Iter_Num, J_SCP,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
% grid on;
% xlabel("SCP Iteration Number",'interpreter','latex','FontSize',15)
% ylabel("Cost Function Value",'interpreter','latex','FontSize',15)
% set(gca, 'FontSize', 15);

% figure('Renderer', 'painters', 'Position', [700 500 550 425]) %#ok<FGREN> 
% boxplot(SCP_time);
% ylabel("$solving~time (s)$",'interpreter','latex','FontSize',12)