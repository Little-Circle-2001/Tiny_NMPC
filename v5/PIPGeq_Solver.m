function [PIPG_x_output, PIPG_u_output, PIPG_lambda_output, IterNum] = PIPGeq_Solver(H, h, G, g, x0, u0, lambda0, Scl_Mat,...
                                                                          PIPG_MaxIter, PIPG_Tol, alpha, beta)
    % 状态量x的维度为(x_dim, N+1)
    % 控制量u的维度为(u_dim, N)
    % 动力学约束拉格朗日乘子λ的维度为(x_dim, N)
    % 
    % 标准Alg.1a (PIPGeq)实现 - 使用向量化更新和比例项
    % 参考: pipgeq_demo.ipynb中的execute_alg1a!函数
    
    StateDim = size(x0,1);
    InputDim = size(u0,1);
    N = size(u0,2);

%     z_opt = quadprog(H,h,[],[],G,g,[],[],[reshape(x0(:,2:end),[],1);reshape(u0,[],1)]);

    %% Compute step sizes using Alg.1a formula if not provided
    % 如果alpha或beta未提供或为NaN，则计算最优步长
    if nargin < 11 || isnan(alpha) || isnan(beta)
        % 计算约束矩阵的谱范数
        sigma = norm(G'*G);  % 或者使用 svd(G'*G) 的最大奇异值
        
        % 计算成本矩阵的最大特征值
        lambda_max = norm(H);  % 或者使用 svd(H) 的最大奇异值
        
        % Alg.1a步长公式: alpha = (-lambda + sqrt(lambda^2 + 4*sigma))/(2*sigma)
        alpha = (-lambda_max + sqrt(lambda_max^2 + 4*sigma))/(2*sigma);
        beta = alpha;  % Alg.1a中 ρ = α
    end

    %% Convert to vectorized form
    % 将x0, u0转换为向量形式 z = [x(2:end); u]
    z = [reshape(x0(:,2:end), [], 1); reshape(u0, [], 1)];
    
    % 将lambda0转换为向量形式 w (只取k=2:N+1的部分)
    w = reshape(lambda0(:,2:end), [], 1);
    
    % 获取矩阵组件用于后续分解
    [Q, S, R, q, r, A, B, gamma] = Matrix_Component(H, h, G, g, StateDim, InputDim, N);

    %% Initialize output
    PIPG_x_output = x0;
    PIPG_u_output = u0;
    PIPG_lambda_output = lambda0;

    %% Begin Alg.1a Iteration (vectorized form)
    p = 1;
    z_prev = z;  % 保存上一次的z值用于收敛判断
    while p <= PIPG_MaxIter
        % Alg.1a核心迭代（向量化形式）
        % 计算约束违反
        constraint_violation = G*z - g;
        
        % 更新z: z = z - alpha*(H*z + h + G'*w + beta*G'*(G*z - g))
        % 注意：beta*G'*(G*z - g) 是比例项（proportional term）
        z = z - alpha*(H*z + h + G'*w + beta*G'*constraint_violation);
        
        % 更新w: w = w + beta*(G*z - g)
        w = w + beta*constraint_violation;
        
        % 将z分解回x和u
        x_new = reshape(z(1:StateDim*N), StateDim, N);
        u_new = reshape(z(StateDim*N+1:end), InputDim, N);
        
        % 组装完整的x (包含初始状态)
        x_full = [x0(:,1), x_new];
        
        % 将w分解回lambda
        lambda_new = reshape(w, StateDim, N);
        lambda_full = [zeros(StateDim, 1), lambda_new];
        
        % 保存输出
        PIPG_x_output = cat(3, PIPG_x_output, x_full);
        PIPG_u_output = cat(3, PIPG_u_output, u_new);
        PIPG_lambda_output = cat(3, PIPG_lambda_output, lambda_full);
        
        % 检查收敛
        delta_z = z - z_prev;
        Dynamic_res_norm(p) = norm(G*z - g);
        
        if Dynamic_res_norm(end) <= PIPG_Tol && norm(delta_z) <= PIPG_Tol
            break;
        end
        
        z_prev = z;  % 更新上一次的z值
        p = p + 1;
    end

%     disp(max(abs(z_opt-z)));
%     if max(abs(z_opt-z)) >= 7e-4
%         1;
%     end
%     if p > 1950
%         figure()
%         hold on;
%         plot(Dynamic_res_norm, 'LineStyle','-', 'LineWidth', 1.5 , 'Color','b');
%         XlabRang = get(gca,'xlim');  % 获取横坐标范围
%         XlabMmin = XlabRang(1);
%         XlabMmax = XlabRang(2);
%         YlabRang = get(gca,'ylim');
%         YlabMmin = YlabRang(1);
%         YlabMmax = PIPG_Tol;
%         FaceAlpha = 0.8;  % 背景
%         Background = fill([XlabMmin XlabMmin XlabMmax XlabMmax],[YlabMmin YlabMmax YlabMmax YlabMmin],'b','FaceColor','#898989','FaceAlpha',FaceAlpha);
%     end
    IterNum = p;
end
