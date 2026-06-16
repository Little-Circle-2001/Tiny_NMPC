close all;
clear; clc;
% 生成可行的参考轨迹（轨迹不是能量最优）
% dbstop if error %在报错前执行断点
% dbstop if warning %在警告前执行断点

%% 初值与参数设置
% MPC 参数设置
T = 20;
dt = 0.2;
tn = T/dt + 1;
N = 7;

% 设置初值
x0 = 7;
y0 = 7;
z0 = 7;
vx0 = 0;
vy0 = 0;
vz0 = 0;
target = [0;0;0];


X = x0; Y = y0; Z = z0;
v = sqrt(vx0^2+vy0^2+vz0^2); a = [];
tic;
for j = 1:tn-N
    yalmip('clear')
    x = sdpvar(1,N+1);
    y = sdpvar(1,N+1);
    z = sdpvar(1,N+1);
    ax = sdpvar(1,N);
    ay = sdpvar(1,N);
    az = sdpvar(1,N);
    vx = sdpvar(1,N+1);
    vy = sdpvar(1,N+1);
    vz = sdpvar(1,N+1);
    
    % Cost function
    J = 0;
    rho = 1e4;
    for i = 1 : N-1
        J = J + 1/2 * ((x(i)*x(i) + y(i)*y(i) +z(i)*z(i)) + (vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)));
        J = J + 1/2 * (ax(i)*ax(i)+ay(i)*ay(i)+az(i)*az(i));
    end
    J = J + 1/2 * (x(N+1)*x(N+1) + y(N+1)*y(N+1) + z(N+1)*z(N+1) + vx(N+1)*vx(N+1) + vy(N+1)*vy(N+1) + vz(N+1)*vz(N+1));
    
    % Constraints
    Constraints = [];
    Constraints = [Constraints; x(1)==x0; y(1)==y0; z(1)==z0];
    Constraints = [Constraints; vx(1)==vx0; vy(1)==vy0; vz(1)==vz0];
    for i = 1 : N
        Constraints = [Constraints; x(i) >= 0; y(i) >= 0; z(i) >=0];
        Constraints = [Constraints; vx(i+1) == vx(i) + ax(i)*dt];
        Constraints = [Constraints; vy(i+1) == vy(i) + ay(i)*dt];
        Constraints = [Constraints; vz(i+1) == vz(i) + az(i)*dt];
        Constraints = [Constraints; x(i+1) == x(i) + vx(i)*dt + 1/2 * ax(i)*dt^2];
        Constraints = [Constraints; y(i+1) == y(i) + vy(i)*dt + 1/2 * ay(i)*dt^2];
        Constraints = [Constraints; z(i+1) == z(i) + vz(i)*dt + 1/2 * az(i)*dt^2];
        Constraints = [Constraints; -1 <= vx(i) <= 1];
        Constraints = [Constraints; -1 <= vy(i) <= 1];
        Constraints = [Constraints; -1 <= vz(i) <= 1];
        Constraints = [Constraints; -1/1.414 <= ax(i) <= 1/1.414];
        Constraints = [Constraints; -1/1.414 <= ay(i) <= 1/1.414];
        Constraints = [Constraints; -1/1.414 <= az(i) <= 1/1.414];
%         Constraints = [Constraints; ((x0-4)^2 + (y0-4)^2 + (z0-4)^2)+2*(x0-4)*(x(i+1)-x0)+2*(y0-4)*(y(i+1)-y0)+2*(z0-4)*(z(i+1)-z0)>=1.5];
        Constraints = [Constraints; ((x(i)-4)^2 + (y(i)-4)^2 + (z(i)-4)^2)>=1.5]; %线性化可能导致无解
    end

    if j == tn - N
        Constraints = [Constraints; [x(N+1);y(N+1);z(N+1)]==target];
        Constraints = [Constraints; [vx(N+1);vy(N+1);vz(N+1)]==[0;0;0]];
    end
    % Solve
    
%     ops = sdpsettings('verbose',1,'solver','mosek');
%     sol = optimize(Constraints,J,ops);
    sol = optimize(Constraints,J);
%     if sol.problem~=0
%         error('problem = 1');
%     end
    
    AX = value(ax); AY = value(ay); AZ = value(az);
    if j ~= tn - N
        ax1 = AX(1);
        ay1 = AY(1);
        az1 = AZ(1);
        vx1 = vx0 + ax1 * dt;
        vy1 = vy0 + ay1 * dt;
        vz1 = vz0 + az1 * dt;
        x1 = x0 + vx0 * dt + 1/2 * ax1 * dt^2;
        y1 = y0 + vy0 * dt + 1/2 * ay1 * dt^2;
        z1 = z0 + vz0 * dt + 1/2 * az1 * dt^2;
        x0 = x1; y0 = y1; z0 = z1; vx0 = vx1; vy0 = vy1; vz0 = vz1; ax0 = ax1; ay0 = ay1; az0 = az1;
    else
        x1 = value(x(2:end));
        y1 = value(y(2:end));
        z1 = value(z(2:end));
        vx1 = value(vx(2:end));
        vy1 = value(vy(2:end));
        vz1 = value(vz(2:end));
        ax1 = AX;
        ay1 = AY;
        az1 = AZ;
    end
    X = [X,x1];
    Y = [Y,y1];
    Z = [Z,z1];
    v = [v,sqrt(vx1.^2+vy1.^2+vz1.^2)];
    a = [a,sqrt(ax1.^2+ay1.^2+az1.^2)];
end
toc;


A = [X;Y;Z];
save("userReference.mat", "A");