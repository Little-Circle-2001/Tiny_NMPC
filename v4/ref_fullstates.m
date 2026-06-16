% 目标状态target_state不能由用户此时随意给定，因为已经通过先前的路径规划确定了
[xr,ur] = estimate_full_states_NED(A, dt, start);
% plot_reconstruction_results(A, states_full, dt);
function [states_full,input_full] = estimate_full_states_NED(trajectory, dt, start_state, params)
    % 在北东地坐标系中估计无人机的完整状态
    %
    % 输入:
    %   trajectory: 3×m矩阵，每列是[x_N; x_E; x_D] (m)
    %   dt: 采样时间 (s)
    %   start_state: 12×1列向量，初始状态
    %   params: 结构体，包含无人机参数
    %
    % 输出:
    %   states_full: 12*m矩阵，完整状态估计
    %   input_full:4*m矩阵，系统输入估计
    
    if nargin < 4
        params = struct();
    end
    if ~isfield(params, 'm')
        params.m = 2.69;
    end
    if ~isfield(params, 'g')
        params.g = 9.81;
    end
    if ~isfield(params, 'Ixx')
        params.Ixx = 0.015;
    end
    if ~isfield(params, 'Iyy')
        params.Iyy = 0.015;
    end
    if ~isfield(params, 'Izz')
        params.Izz = 0.0245;
    end
    if ~isfield(params, 'smooth_window')
        params.smooth_window = 5; % 平滑窗口大小(此处用不上，后续省略)
    end
    m = params.m;
    g = params.g;
    Ixx = params.Ixx;
    Iyy = params.Iyy;
    Izz = params.Izz;
    smooth_win = params.smooth_window;
    % 获取轨迹维度
    [dim, n_points] = size(trajectory);
    if dim ~= 3
        error('轨迹必须是3×m矩阵，每列是[x_N; x_E; x_D]');
    end
    states_full = zeros(12, n_points);
    input_full = zeros(4,n_points-1);
    % 1. 设置位置
    states_full(1:3, :) = trajectory;
    v_N = zeros(1, n_points);
    v_E = zeros(1, n_points);
    v_D = zeros(1, n_points);
    for i = 1:n_points-1
        v_N(i) = (trajectory(1, i+1) - trajectory(1, i)) / (dt);
        v_E(i) = (trajectory(2, i+1) - trajectory(2, i)) / (dt);
        v_D(i) = (trajectory(3, i+1) - trajectory(3, i)) / (dt);
    end
    
    v_N(end) = (trajectory(1, end) - trajectory(1, end-1)) / dt;
    v_E(end) = (trajectory(2, end) - trajectory(2, end-1)) / dt;
    v_D(end) = (trajectory(3, end) - trajectory(3, end-1)) / dt;
    
%     % 平滑速度数据
%     if smooth_win > 1
%         v_N = movmean(v_N, smooth_win);
%         v_E = movmean(v_E, smooth_win);
%         v_D = movmean(v_D, smooth_win);
%     end
    
    states_full(4:6, :) = [v_N; v_E; v_D];
    
    % 3. 计算偏航角（跟随速度方向）
    psi = zeros(1, n_points);
    for i = 1:n_points
        if norm([v_N(i), v_E(i)]) > 0.01
            psi(i) = atan2(v_E(i), v_N(i)); 
        else
            if i > 1
                psi(i) = psi(i-1);
            else
                psi(i) = start_state(9);
            end
        end
    end
    states_full(9, :) = psi;
    
    % 4. 数值微分计算加速度
    a_N = zeros(1, n_points);
    a_E = zeros(1, n_points);
    a_D = zeros(1, n_points);
    for i = 1:n_points-1
        a_N(i) = (v_N(i+1) - v_N(i)) / (dt);
        a_E(i) = (v_E(i+1) - v_E(i)) / (dt);
        a_D(i) = (v_D(i+1) - v_D(i)) / (dt);
    end
    a_N(end) = (v_N(end) - v_N(end-1)) / dt;
    a_E(end) = (v_E(end) - v_E(end-1)) / dt;
    a_D(end) = (v_D(end) - v_D(end-1)) / dt;
%     % 平滑加速度
%     if smooth_win > 1
%         a_N = movmean(a_N, smooth_win);
%         a_E = movmean(a_E, smooth_win);
%         a_D = movmean(a_D, smooth_win);
%     end
    % 5. 计算z_B和总升力T
    % z_B = - (a - g) / ||a - g||
    % T = m * ||a - g||
    z_B = zeros(3, n_points);
    T = zeros(1, n_points);
    for i = 1:n_points
        a_vec = [a_N(i); a_E(i); a_D(i)];
        a_minus_g = a_vec - [0; 0; g];
        norm_a = norm(a_minus_g);
        if norm_a > 1e-6
            T(i) = m * norm_a;

            z_B(:, i) = -a_minus_g / norm_a;
        else
            T(i) = m * g;
            z_B(:, i) = [0; 0; 1];
        end
    end
    % 6. 计算滚转角φ和俯仰角θ
    phi = zeros(1, n_points);
    theta = zeros(1, n_points);
    for i = 1:n_points
        z1 = z_B(1, i);
        z2 = z_B(2, i);
        z3 = z_B(3, i);
        psi_i = psi(i);
        A = z1 * cos(psi_i) + z2 * sin(psi_i);
        B = -z1 * sin(psi_i) + z2 * cos(psi_i);
        C = z3;
        if abs(B) < 0.9999
            phi(i) = -asin(B);
        else
            phi(i) = -sign(B) * pi/2;
        end
        % θ = atan2(A, C)
        if abs(A) > 1e-6 || abs(C) > 1e-6
            theta(i) = atan2(A, C);
        else
            theta(i) = 0;
        end
        phi(i) = atan2(sin(phi(i)), cos(phi(i)));
        theta(i) = atan2(sin(theta(i)), cos(theta(i)));
    end
    states_full(7, :) = phi;    % 滚转角
    states_full(8, :) = theta;  % 俯仰角
    
    % 7. 计算角速度p, q, r
    % 首先计算欧拉角的时间导数
    phi_dot = zeros(1, n_points);
    theta_dot = zeros(1, n_points);
    psi_dot = zeros(1, n_points);
    
    for i = 1:n_points-1
        phi_dot(i) = (phi(i+1) - phi(i)) / (dt);
        theta_dot(i) = (theta(i+1) - theta(i)) / (dt);
        psi_dot(i) = (psi(i+1) - psi(i)) / (dt);
    end
    phi_dot(end) = (phi(end) - phi(end-1)) / dt;
    theta_dot(end) = (theta(end) - theta(end-1)) / dt;
    psi_dot(end) = (psi(end) - psi(end-1)) / dt;
%     % 平滑角速度
%     if smooth_win > 1
%         phi_dot = movmean(phi_dot, smooth_win);
%         theta_dot = movmean(theta_dot, smooth_win);
%         psi_dot = movmean(psi_dot, smooth_win);
%     end
    % 将欧拉角速度转换为机体角速度
    p = zeros(1, n_points);
    q = zeros(1, n_points);
    r = zeros(1, n_points);
    
    for i = 1:n_points
        phi_i = phi(i);
        theta_i = theta(i);
        W_inv = [1, 0, -sin(theta_i);
                 0, cos(phi_i), sin(phi_i)*cos(theta_i);
                 0, -sin(phi_i), cos(phi_i)*cos(theta_i)];
        omega = W_inv * [phi_dot(i); theta_dot(i); psi_dot(i)];
        p(i) = omega(1);
        q(i) = omega(2);
        r(i) = omega(3);
    end
    states_full(10, :) = p; 
    states_full(11, :) = q; 
    states_full(12, :) = r; 
    % 8. 用初始状态修正第一个点的状态
    states_full(:, 1) = start_state;
    % 9. 计算角加速度p_dot,q_dot,r_dot
    p_dot = zeros(1, n_points);
    q_dot = zeros(1, n_points);
    r_dot = zeros(1, n_points);
    for i = 1:n_points-1
        p_dot(i) = (p(i+1) - p(i)) / (dt);
        q_dot(i) = (q(i+1) - q(i)) / (dt);
        r_dot(i) = (r(i+1) - r(i)) / (dt);
    end
    p_dot(end) = (p(end) - p(end-1)) / dt;
    q_dot(end) = (q(end) - q(end-1)) / dt;
    r_dot(end) = (r(end) - r(end-1)) / dt;
    % 10.计算参考力矩
    input_full(1,:) = T(1:end-1);
    for i = 1:n_points-1
        input_full(2,i) = p_dot(i)*Ixx-(Iyy-Izz)*q(i)*r(i);
        input_full(3,i) = q_dot(i)*Iyy-(Izz-Ixx)*p(i)*r(i);
        input_full(4,i) = r_dot(i)*Izz-(Ixx-Iyy)*p(i)*q(i);
    end
end
function plot_reconstruction_results(trajectory, states_full, dt)
    % 可视化重构结果
    
    n_points = size(trajectory, 2);
    t = (0:n_points-1) * dt;
    
    figure('Position', [100, 100, 1200, 800]);
    
    % 1. 3D轨迹
    subplot(2, 3, 1);
    plot3(trajectory(1, :), trajectory(2, :), trajectory(3, :));
    hold on;
    plot3(trajectory(1, 1), trajectory(2, 1), trajectory(3, 1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot3(trajectory(1, end), trajectory(2, end), trajectory(3, end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('北向 (m)'); ylabel('东向 (m)'); zlabel('高度 (m)');
    title('3D轨迹（高度向上为正）');
    grid on; view(45, 30);
    axis equal;
    
    % 2. 位置
    subplot(2, 3, 2);
    plot(t, trajectory(1, :), 'b-', 'LineWidth', 2); hold on;
    plot(t, trajectory(2, :), 'r-', 'LineWidth', 2);
    plot(t, trajectory(3, :), 'g-', 'LineWidth', 2);
    xlabel('时间 (s)'); ylabel('位置 (m)');
    title('位置分量');
    legend('北向', '东向', '高度');
    grid on;
    
    % 3. 速度
    subplot(2, 3, 3);
    plot(t, states_full(4, :), 'b-', 'LineWidth', 2); hold on;
    plot(t, states_full(5, :), 'r-', 'LineWidth', 2);
    plot(t, states_full(6, :), 'g-', 'LineWidth', 2);
    xlabel('时间 (s)'); ylabel('速度 (m/s)');
    title('速度分量');
    legend('v_N', 'v_E', 'v_D');
    grid on;
    
    % 4. 姿态角
    subplot(2, 3, 4);
    plot(t, rad2deg(states_full(7, :)), 'b-', 'LineWidth', 2); hold on;
    plot(t, rad2deg(states_full(8, :)), 'r-', 'LineWidth', 2);
    plot(t, rad2deg(states_full(9, :)), 'g-', 'LineWidth', 2);
    xlabel('时间 (s)'); ylabel('角度 (°)');
    title('姿态角');
    legend('滚转φ', '俯仰θ', '偏航ψ');
    grid on;
    
    % 5. 角速度
    subplot(2, 3, 5);
    plot(t, rad2deg(states_full(10, :)), 'b-', 'LineWidth', 2); hold on;
    plot(t, rad2deg(states_full(11, :)), 'r-', 'LineWidth', 2);
    plot(t, rad2deg(states_full(12, :)), 'g-', 'LineWidth', 2);
    xlabel('时间 (s)'); ylabel('角速度 (°/s)');
    title('机体角速度');
    legend('p', 'q', 'r');
    grid on;
    
    % 6. 速度方向
    subplot(2, 3, 6);
    speed = sqrt(states_full(4, :).^2 + states_full(5, :).^2 + states_full(6, :).^2);
    plot(t, speed, 'b-', 'LineWidth', 2);
    xlabel('时间 (s)'); ylabel('速度大小 (m/s)');
    title('总速度');
    grid on;
    
    sgtitle('状态重构结果');
end