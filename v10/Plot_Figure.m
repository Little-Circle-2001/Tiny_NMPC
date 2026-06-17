figure(1)
hold on;
% plot obstacle
th = linspace(0,2*pi*100);
x = cos(th) ; y = sin(th) ;
obstacle_Xedge = []; obstacle_Yedge = [];

for i = 1:1
    obstacle_Xedge(:,i) = (4 + 1*x)'; %#ok<SAGROW> 
    obstacle_Yedge(:,i) = (4 + 1*y)'; %#ok<SAGROW> 
%     l4(i) = plot(obstacle_Xedge(:,i), obstacle_Yedge(:,i), 'Color', [0.8500, 0.3250, 0.0980],...
%         'LineWidth', 2);
%     l5(i) = plot(pos(1,i), pos(2,i), 'Color', [0.8500, 0.3250, 0.0980],...
%         'MarkerSize', 5, 'LineWidth', 2);
end
Obstacle = fill(obstacle_Xedge,obstacle_Yedge,'white','EdgeColor',[0.8500, 0.3250, 0.0980],'FaceAlpha',0.9);
plot(state_SCP_opt(1,:,end),state_SCP_opt(2,:,end),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
grid on;
xlabel("$x_1 (m)$",'interpreter','latex','FontSize',12)
ylabel("$x_2 (m)$",'interpreter','latex','FontSize',12)
set(gca, 'FontSize', 12);

figure(2)
plot(input_SCP_opt(1,:,end)*58.3,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
hold on;
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$\delta (deg)$",'interpreter','latex','FontSize',12)
legend("$\delta$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);

figure(3)
hold on;
plot(input_SCP_opt(2,:,end),'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','b')
grid on;
xlabel("$t (s)$",'interpreter','latex','FontSize',12)
ylabel("$a (m \cdot s^{-2})$",'interpreter','latex','FontSize',12)
legend("$a$",'interpreter','latex','FontSize',12);
set(gca, 'FontSize', 12);

figure('Renderer', 'painters', 'Position', [200 100 1200 600]) %#ok<FGREN> 
plot(0:SCP_Iter_Num, J_SCP,'Marker','o','LineStyle','-', 'LineWidth', 1.5 , 'Color','r')
grid on;
xlabel("SCP Iteration Number",'interpreter','latex','FontSize',15)
ylabel("Cost Function Value",'interpreter','latex','FontSize',15)
set(gca, 'FontSize', 15);