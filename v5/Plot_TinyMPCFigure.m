close all
figure('Renderer', 'painters', 'Position', [600 100 800 700]) %#ok<FGREN>

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

h1 = plot3(xlog(1,:),xlog(2,:),xlog(3,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r');
h2 = plot3(xr(1,:),xr(2,:),xr(3,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b');
grid on;
xlabel("$x_1 (m)$",'interpreter','latex','FontSize',12)
ylabel("$x_2 (m)$",'interpreter','latex','FontSize',12)
legend([h1,h2],"Actual","Reference",'interpreter','latex','FontSize',12)
set(gca, 'FontSize', 12);

figure('Renderer', 'painters', 'Position', [700 500 1100 425]) %#ok<FGREN>
subplot(1,2,1)
plot(xlog(4,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
hold on;
plot(xlog(5,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','g')
plot(xlog(6,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','b')
% plot(xr(4,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b');
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$v (m/s)$",'interpreter','latex','FontSize',12)
legend("$v_x$","$v_y$","$v_z$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);

subplot(1,2,2)
plot(xr(4,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','r');
hold on;
plot(xr(5,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','g');
plot(xr(6,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b');
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$v_r (m/s)$",'interpreter','latex','FontSize',12)
legend("$v_{rx}$","$v_{ry}$","$v_{rz}$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);
figure('Renderer', 'painters', 'Position', [700 500 1100 425]) %#ok<FGREN> 
subplot(1,2,1)
hold on;
plot(ulog(1,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
plot(ulog(2,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','g')
plot(ulog(3,:),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','b')
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$u (m/s^2)$",'interpreter','latex','FontSize',12)
legend("$u_x$","$u_y$","$u_z$",'interpreter','latex','FontSize',12);

subplot(1,2,2)
hold on;
plot(ur(1,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','r')
plot(ur(2,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','g')
plot(ur(3,:),'LineStyle','--', 'LineWidth', 1.5 , 'Color','b')
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$u_r (m/s^2)$",'interpreter','latex','FontSize',12)
legend("$u_{rx}$","$u_{ry}$","$u_{rz}$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);

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

figure('Renderer', 'painters', 'Position', [700 500 550 425]) %#ok<FGREN> 
boxplot(SCP_time);
ylabel("$solving~time (s)$",'interpreter','latex','FontSize',12)