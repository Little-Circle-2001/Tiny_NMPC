function [PIPG_x_output, PIPG_u_output, PIPG_lambda_output, IterNum] = PIPG_Solver(H, h, G, g, x0, u0, lambda0, Scl_Mat, ...
                                                                          PIPG_MaxIter, PIPG_Tol, alpha, beta)
    % 状态量x的维度为(x_dim, N+1)
    % 控制量u的维度为(u_dim, N)
    % 动力学约束拉格朗日乘子λ的维度为(x_dim, N)
    
    StateDim = size(x0,1);
    InputDim = size(u0,1);
    N = size(u0,2);
    z0 = [reshape(x0(:,2:end), [], 1);...
          reshape(u0, [], 1)];
    lambda_vec0 = reshape(lambda0(:,2:end), [], 1);

    z_opt = quadprog(H,h,[],[],G,g,[],[],[reshape(x0(:,2:end),[],1);reshape(u0,[],1)]);

    H0 = eye(size(z0,1));

    %% Initialize
    p = 1;

    PIPG_x_output = x0;
    PIPG_u_output = u0;
    PIPG_lambda_output = lambda0;


    %% Begin PIPG Iteration
%     tic;
    while p <= PIPG_MaxIter

        xi = lambda_vec0 + beta * (G * z0 - g);

        z1 = z0 - alpha * (H * z0 + h + G' * xi);

%         grad_f0 = H * z0 + h + G' * xi;
%         p0 = -H0 * grad_f0;
%         z1 = z0 + alpha * p0;
%         grad_f1 = H * z1 + h + G' * xi;
%         s = z1 - z0;
%         y = grad_f1 - grad_f0;
%         rho = 1/(y'*s);
%         H1 = (eye(length(s))-rho*s*y') * H0 * (eye(length(s))-rho*y*s') + rho*s*s';
%         H0 = H1;
%         z0 = z1;
%         grad_f0 = grad_f1;

        Act_z = Scl_Mat * z1;

%         lambda_vec1 = lambda_vec0 + beta * (G * z1 - g);
        lambda_vec1 = xi + beta * G * (z1 - z0);

        x1 = [x0(:,1), reshape(Act_z(1:StateDim*N),StateDim,N)];
        u1 = [reshape(Act_z(StateDim*N+1:end),InputDim,N)];
        lambda1 = [lambda0(:,1), reshape(lambda_vec1,StateDim,N)];

        PIPG_x_output = cat(3, PIPG_x_output, x1);
        PIPG_u_output = cat(3, PIPG_u_output, u1);
        PIPG_lambda_output = cat(3, PIPG_lambda_output, lambda1);

        delta_z = z1 - z0;
        Dynamic_res_norm(p) = norm(G*z1-g);

        z0 = z1;
        lambda_vec0 = lambda_vec1;
        
%         if Dynamic_res_norm(end) <= PIPG_Tol && norm(delta_z) <= PIPG_Tol
%             break;
%         end
        if Dynamic_res_norm(end) <= PIPG_Tol
            break;
        end

        if p == 1000
            1;
        end
        p = p + 1;
    end

%     disp(G*z1-g)
%     disp(max(abs(z_opt-z1)));
%     if max(abs(z_opt-z)) >= 7e-4
%         1;
%     end


%     if p > 2000
%         figure()
%         hold on;
%         plot(Dynamic_res_norm, 'LineStyle','-', 'LineWidth', 1.5 , 'Color','b');
%         XlabRang = get(gca,'xlim');  % 获取横坐标范围
%         XlabMmin = XlabRang(1);
%         XlabMmax = XlabRang(2);
%         YlabRang = get(gca,'ylim');
%         YlabMmin = 0;
%         YlabMmax = PIPG_Tol;
%         FaceAlpha = 0.8;  % 背景
%         Background = fill([XlabMmin XlabMmin XlabMmax XlabMmax],[YlabMmin YlabMmax YlabMmax YlabMmin],'b','FaceColor','#898989','FaceAlpha',FaceAlpha);
%     end

    IterNum = p;
end