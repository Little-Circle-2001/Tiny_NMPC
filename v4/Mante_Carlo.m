close all;
clear; clc;
global tn; %#ok<*GVMIS> 

%% 初值与参数设置
% MPC 参数设置
N = 10;                      % 预测时域
T = 600;
dt = 6;
tn = T/dt + 1;

mu = 3.986e5;
a = 6871.393;
n = sqrt(mu/a^3);

% 参考轨迹 设置初值
%% 设置初值
random = normrnd(0,0.025,128,3);
x_mante_log = [];
Monte_Time =  [];
Monte_avgTime = [];
Mosek_Monte_Time =  [];
Mosek_Monte_avgTime = [];
OSQP_Mante_Time =  [];
OSQP_Monte_avgTime = [];
for Monte_num = 1 : size(random,1)
    x0 = -0.500 + random(Monte_num,1);
    y0 = 1.000 + random(Monte_num,2);
    z0 = 0.500 + random(Monte_num,3);
    vx0 = 0;
    vy0 = 0;
    vz0 = 0;
    % start0 = [3;0;-pi;0.3];
    target = [0;0;0];
    
    %% Obstacle
    Obs_Pos = [-0.2 -0.4; ...
                0.4  0.7; ...
                0.15  0.35];
    Obs_r = [0.1; 0.1];

    
    X = x0; Y = y0; Z = z0;
    VX = vx0; VY = vy0; VZ = vz0;
    AX = []; AY = []; AZ = [];
    
    %%
%     tic;

    for j = 1:tn-N
        yalmip('clear')
        x = sdpvar(1,N+1);
        y = sdpvar(1,N+1);
        z = sdpvar(1,N+1);
        ux = sdpvar(1,N);
        uy = sdpvar(1,N);
        uz = sdpvar(1,N);
        vx = sdpvar(1,N+1);
        vy = sdpvar(1,N+1);
        vz = sdpvar(1,N+1);
        
        % Cost function
        J = 0;
        for i = 1 : N
            J = J + 1/2 * 1e-10 * (x(i)*x(i) + y(i)*y(i) + z(i)*z(i));
            J = J + 1/2 * 50 * (vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i));
            J = J + 1/2 * 50 * (uy(i)*uy(i) + uy(i)*uy(i) + uz(i)*uz(i));
        end
        J = J + 1/2 * (x(N+1)*x(N+1) + y(N+1)*y(N+1) + z(N+1)*z(N+1) + vx(N+1)*vx(N+1) + vy(N+1)*vy(N+1) + vz(N+1)*vz(N+1));
        
        % Constraints
        Constraints = [];
        Constraints = [Constraints; [x(1);y(1);z(1)]==[x0;y0;z0]];
        Constraints = [Constraints; vx(1)==vx0; vy(1)==vy0; vz(1)==vz0];
        for i = 1 : N
            Constraints = [Constraints; x(i) <= 0; y(i) >= 0; z(i) >= 0];
            Constraints = [Constraints; vx(i+1) == vx(i) + ux(i)*dt];
            Constraints = [Constraints; vy(i+1) == vy(i) + uy(i)*dt];
            Constraints = [Constraints; vz(i+1) == vz(i) + uz(i)*dt];
            Constraints = [Constraints; x(i+1) == x(i) + vx(i)*dt + 1/2 * ux(i)*dt^2];
            Constraints = [Constraints; y(i+1) == y(i) + vy(i)*dt + 1/2 * uy(i)*dt^2];
            Constraints = [Constraints; z(i+1) == z(i) + vz(i)*dt + 1/2 * uz(i)*dt^2];
            Constraints = [Constraints; -0.004 <= vx(i) <= 0.004];
            Constraints = [Constraints; -0.004 <= vy(i) <= 0.004];
            Constraints = [Constraints; -0.004 <= vz(i) <= 0.004];
            Constraints = [Constraints; -0.00008 <= ux(i) <= 0.00008];
            Constraints = [Constraints; -0.00008 <= uy(i) <= 0.00008];
            Constraints = [Constraints; -0.00008 <= uz(i) <= 0.00008];
            for k = 1 : length(Obs_r)
                if (x0-Obs_Pos(1,k))^2+(y0-Obs_Pos(2,k))^2+(z0-Obs_Pos(3,k))^2 <= 0.3^2
                    Constraints = [Constraints; (x0-Obs_Pos(1,k))^2+(y0-Obs_Pos(2,k))^2+(z0-Obs_Pos(3,k))^2+ ...
                                                2*(x0-Obs_Pos(1,k))*(x(i+1)-x0)+ ...
                                                2*(y0-Obs_Pos(2,k))*(y(i+1)-y0)+ ...
                                                2*(z0-Obs_Pos(3,k))*(z(i+1)-z0) >= Obs_r(k)^2];
                end
            end
        end
    
        if j == tn - N
            Constraints = [Constraints; [x(N+1);y(N+1);z(N+1)]==target];
            Constraints = [Constraints; [vx(N+1);vy(N+1);vz(N+1)]==[0;0;0]];
        end
        % Solve
        
        ops = sdpsettings('verbose',0,'solver','mosek');
        sol = optimize(Constraints,J,ops);
        
        UX = value(ux); UY = value(uy); UZ = value(uz);
        if j ~= tn - N
            ux1 = UX(1);
            uy1 = UY(1);
            uz1 = UZ(1);
            vx1 = vx0 + ux1 * dt;
            vy1 = vy0 + uy1 * dt;
            vz1 = vz0 + uz1 * dt;
            x1 = x0 + vx0 * dt + 1/2 * ux1 * dt^2;
            y1 = y0 + vy0 * dt + 1/2 * uy1 * dt^2;
            z1 = z0 + vz0 * dt + 1/2 * uz1 * dt^2;
            x0 = x1; y0 = y1; z0 = z1;
            vx0 = vx1; vy0 = vy1; vz0 = vz1;
            ux0 = ux1; uy0 = uy1; uz0 = uz1;
        else
            x1 = value(x(2:end));
            y1 = value(y(2:end));
            z1 = value(z(2:end));
            vx1 = value(vx(2:end));
            vy1 = value(vy(2:end));
            vz1 = value(vz(2:end));
            ux1 = UX;
            uy1 = UY;
            uz1 = UZ;
        end
        if sol.problem == 0
            X = [X,x1];
            Y = [Y,y1];
            Z = [Z,z1];
            VX = [VX,vx1];
            VY = [VY,vy1];
            VZ = [VZ,vz1];
            AX = [AX,ux1];
            AY = [AY,uy1];
            AZ = [AZ,uz1];
        else
            break;
        end
    
    end
%     toc;
    
    
    A = [X;Y;Z;VX;VY;VZ];
    B = [AX;AY;AZ];
    
    %% 在线求解
    start = [-0.5 + random(Monte_num,1);1 + random(Monte_num,2);0.5 + random(Monte_num,3);0;0;0]; % -500;1000;500;0;0;0
    target = [0;0;0;0;0;0];
    
    xr = A;
    ur = B;
     
    % SCP 求解器参数设置
    SCP_MaxIter = 20;           % SCP最大迭代次数
    SCP_Tol = 1e-5;             % SCP迭代容许收敛误差
    
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
%         fprintf('ADMM运算时间为%.6f秒\n',SCP_time(i))
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
    x_mante_log = cat(3, x_mante_log, xlog);
    Monte_Time(Monte_num) =  sum(SCP_time);
    Monte_avgTime(Monte_num) = mean(SCP_time);
    fprintf('第%d次蒙特卡洛模拟完成，用时%.6f秒, 平均单次SCP用时%.6f秒\n', Monte_num, Monte_Time(Monte_num), Monte_avgTime(Monte_num))

    %% Mosek求解器

    state0 = xr(:,1:N+1);
    input0 = ur(:,1:N);
    Mosek_xlog = start;
    Mosek_ulog = [];
    
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
    
        Mosek_ulog = [Mosek_ulog, u]; %#ok<*AGROW> 
        Mosek_xlog = [Mosek_xlog, x];
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
    
        Mosek_ulog = [Mosek_ulog, u]; %#ok<*AGROW> 
        Mosek_xlog = [Mosek_xlog, x];
    end
    Mosek_Monte_Time(Monte_num) =  sum(Mosek_Time);
    Mosek_Monte_avgTime(Monte_num) = mean(Mosek_Time);
    fprintf('Mosek第%d次蒙特卡洛模拟完成，用时%.6f秒, 平均单次SCP用时%.6f秒\n', Monte_num, Mosek_Monte_Time(Monte_num), Mosek_Monte_avgTime(Monte_num))
    %% OSQP求解器
    
    state0 = xr(:,1:N+1);
    input0 = ur(:,1:N);
    OSQP_xlog = start;
    OSQP_ulog = [];
    
    for j = tn-N+1 : tn-1
        u = input_SCP_opt(:,j-tn+N+1,end);
        x = RK_dyn(x,u,dt,n);
    
        OSQP_ulog = [OSQP_ulog, u]; %#ok<*AGROW> 
        OSQP_xlog = [OSQP_xlog, x];
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
    
        OSQP_ulog = [OSQP_ulog, u]; %#ok<*AGROW> 
        OSQP_xlog = [OSQP_xlog, x];
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
    
        OSQP_ulog = [OSQP_ulog, u]; %#ok<*AGROW> 
        OSQP_xlog = [OSQP_xlog, x];
    end
    OSQP_Monte_Time(Monte_num) =  sum(OSQP_Time);
    OSQP_Monte_avgTime(Monte_num) = mean(OSQP_Time);
    fprintf('OSQP第%d次蒙特卡洛模拟完成，用时%.6f秒, 平均单次SCP用时%.6f秒\n', Monte_num, OSQP_Monte_Time(Monte_num), OSQP_Monte_avgTime(Monte_num))
end
save("x_mante_log.mat", "x_mante_log")
save("Monte_avgTime.mat", "Monte_avgTime")
save("Mosek_Monta_avgTime.mat", "Mosek_Monte_avgTime")
save("OSQP_Monta_avgTime.mat", "OSQP_Monte_avgTime")

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
for i = 1 : size(random,1)
    h1 = plot3(1000 * x_mante_log(1,:,i),1000 * x_mante_log(2,:,i),1000 * x_mante_log(3,:,i),'LineStyle','-', 'LineWidth', 1 , 'Color', rand(1, 3));
    grid on;
    xlabel("$x_1 (m)$",'interpreter','latex','FontSize',12)
    ylabel("$x_2 (m)$",'interpreter','latex','FontSize',12)
    set(gca, 'FontSize', 12);
end
view(48,34)
save2pdf('figures/ManteCarloTraj.pdf')

rng default;
figure('Renderer', 'painters', 'Position', [700 500 550 425]) %#ok<FGREN> 
% group = categorical({'Proposed'; 'Mosek'; 'OSQP'}); % 条件标签
hold on

boxchart(ones(size(Monte_avgTime)), Monte_avgTime, ...
    'BoxFaceColor',[0 0.5 0.5])

boxchart(2*ones(size(Mosek_Monte_avgTime)), Mosek_Monte_avgTime, ...
    'BoxFaceColor',[0.5 0 0.5])

boxchart(3*ones(size(OSQP_Monte_avgTime)), OSQP_Monte_avgTime, ...
    'BoxFaceColor',[0.5 0 0])
xticks([1 2 3]);
xticklabels({'Proposed','Mosek','OSQP'})
set(gca,'FontSize',12)
xlabel('$Solvers$','interpreter','latex','FontSize',12)
ylabel("$Solving~Time (s)$",'interpreter','latex','FontSize',12)
grid on;
save2pdf('figures/Monte_Carlo_ComparedTime.pdf')