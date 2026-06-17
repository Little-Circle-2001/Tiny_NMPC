close all;
clear; clc;

%% 初值与参数设置
% MPC 参数设置
T = 600;
dt = 6;
tn = T/dt + 1;
N = 10;                      % 预测时域

mu = 3.986e5;
a = 6871.393;
n = sqrt(mu/a^3);

%% Obstacle
Obs_Pos = [-0.2 -0.4; ...
            0.4  0.7; ...
            0.15 0.35];
Obs_r = [0.100; 0.100];

%% 设置初值
x0 = -0.500;
y0 = 1.000;
z0 = 0.500;
vx0 = 0;
vy0 = 0;
vz0 = 0;
% start0 = [3;0;-pi;0.3];
target = [0;0;0];


X = x0; Y = y0; Z = z0;
VX = vx0; VY = vy0; VZ = vz0;
AX = []; AY = []; AZ = [];

%%
tic;
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
toc;


A = [X;Y;Z;VX;VY;VZ];
B = [AX;AY;AZ];
save("ReferenceState.mat", "A")
save("ReferenceInput.mat", "B")
% tic
% x = quadprog(H,h,C,-c,G,g,[],[],[],[]);
% fprintf('自带的求解器运算时间为%.6f秒\n',toc);

figure(1)
% plot obstacle
for k = 1 : length(Obs_r)
    [obs_x, obs_y, obs_z] = sphere(50); % 50是网格的分辨率
    obs_x = Obs_r(k) * obs_x + Obs_Pos(1,k);
    obs_y = Obs_r(k) * obs_y + Obs_Pos(2,k);
    obs_z = Obs_r(k) * obs_z + Obs_Pos(3,k);
    
    surf(obs_x, obs_y, obs_z); % 使用'none'去除边缘线
    colormap(hot); % 使用hot颜色映射
    
    axis equal; % 保持轴比例相等，使球体看起来是圆形的
    
    hold on;
end

plot3(X,Y,Z,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
grid on;
xlabel("$x_1 (m)$",'interpreter','latex','FontSize',12)
ylabel("$x_2 (m)$",'interpreter','latex','FontSize',12)
zlabel("$x_3 (m)$",'interpreter','latex','FontSize',12)
set(gca, 'FontSize', 12);

figure(2)
hold on;
plot(VX,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
plot(VY,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','g')
plot(VZ,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','b')
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$v (m/s)$",'interpreter','latex','FontSize',12)
legend("$v_x$","$v_y$","$v_z$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);

figure(3)
hold on;
plot(AX,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
plot(AY,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','g')
plot(AZ,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','b')
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$u (m/s^2)$",'interpreter','latex','FontSize',12)
legend("$u_x$","$u_y$","$u_z$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);
