function [PIPG_x_output, PIPG_u_output, PIPG_lambda_output] = PIPG_Solver_WithProportional(H, h, G, g, x0, u0, lambda0, ...
                                                                          PIPG_MaxIter, PIPG_Tol, alpha, beta)
    % 状态量x的维度为(x_dim, N+1)
    % 控制量u的维度为(u_dim, N)
    % 动力学约束拉格朗日乘子λ的维度为(x_dim, N)
    % 
    % 此版本在x和u更新中添加了比例项（proportional term），用于验证比例项对收敛速度的影响
    
    StateDim = size(x0,1);
    InputDim = size(u0,1);
    N = size(u0,2);

    z_opt = quadprog(H,h,[],[],G,g,[],[],[reshape(x0(:,2:end),[],1);reshape(u0,[],1)]);

    %% Initialize
    [Q, S, R, q, r, A, B, gamma] = Matrix_Component(H, h, G, g, StateDim, InputDim, N);
    p = 1;
    lambda = zeros(StateDim, N+1, PIPG_MaxIter+1);
    x = zeros(StateDim, N+1, PIPG_MaxIter+1);
    u = zeros(InputDim, N, PIPG_MaxIter+1);
    xi = zeros(StateDim, N+1, PIPG_MaxIter);

    x(:,:,1) = x0;
    u(:,:,1) = u0;
    lambda(:,:,1) = lambda0;

    PIPG_x_output = x0;
    PIPG_u_output = u0;
    PIPG_lambda_output = lambda0;

    %% Begin PIPG Iteration
    while p <= PIPG_MaxIter
        for k = 2 : N+1
            if k == 2
                xi(:,k,p) = lambda(:,k,p) + beta * (...
                            x(:,k,p) - B(:,:,k-1)*u(:,k-1,p) - gamma(:,k)...
                            );
            else
                xi(:,k,p) = lambda(:,k,p) + beta * (...
                            x(:,k,p) - B(:,:,k-1)*u(:,k-1,p) - A(:,:,k-1)*x(:,k-1,p) - gamma(:,k)...
                            );
            end
        end

        x(:,1,p+1) = x0(:,1);
        for k = 2 : N+1
            % 计算约束违反项（用于比例项）
            if k == 2
                constraint_violation = x(:,k,p) - B(:,:,k-1)*u(:,k-1,p) - gamma(:,k);
            else
                constraint_violation = x(:,k,p) - B(:,:,k-1)*u(:,k-1,p) - A(:,:,k-1)*x(:,k-1,p) - gamma(:,k);
            end
            
            if k <= N                
                % 添加比例项: + beta * constraint_violation
                x(:,k,p+1) = x(:,k,p) - alpha * ( ...
                             Q(:,:,k)*x(:,k,p) + S(:,:,k)*u(:,k,p) + q(:,k) + xi(:,k,p) - A(:,:,k)'*xi(:,k+1,p) ...
                             + beta * constraint_violation ...  % 比例项
                             );
            else
                x(:,k,p+1) = x(:,k,p) - alpha * ( ...
                             Q(:,:,k)*x(:,k,p) + q(:,k) + xi(:,k,p) ...
                             + beta * constraint_violation ...  % 比例项
                             );
            end
        end
        for k = 2 : N+1
            % 计算约束违反项（用于比例项）
            if k == 2
                constraint_violation = x(:,k,p) - B(:,:,k-1)*u(:,k-1,p) - gamma(:,k);
            else
                constraint_violation = x(:,k,p) - B(:,:,k-1)*u(:,k-1,p) - A(:,:,k-1)*x(:,k-1,p) - gamma(:,k);
            end
            
            % 添加比例项: + beta * constraint_violation
            u(:,k-1,p+1) = u(:,k-1,p) - alpha * (...
                           S(:,:,k-1)'*x(:,k-1,p) + R(:,:,k-1)*u(:,k-1,p) + r(:,k-1) - B(:,:,k-1)'*xi(:,k,p)...
                           + beta * constraint_violation ...  % 比例项
                           );
        end
        for k = 2 : N+1
            lambda(:,k,p+1) = lambda(:,k,p) + beta * (...
                              x(:,k,p+1) - B(:,:,k-1)*u(:,k-1,p+1) - A(:,:,k-1)*x(:,k-1,p+1) - gamma(:,k)...
                              );
        end
        PIPG_x_output = cat(3, PIPG_x_output, x(:,:,p+1));
        PIPG_u_output = cat(3, PIPG_u_output, u(:,:,p+1));
        PIPG_lambda_output = cat(3, PIPG_lambda_output, lambda(:,:,p+1));

        z = [reshape(x(:,2:end,p+1), [], 1);...
             reshape(u(:,:,p+1), [], 1)];
        delta_x = x(:,2:end,p+1) - x(:,2:end,p);
        delta_u = u(:,:,p+1) - u(:,:,p);
        Dynamic_res_norm = norm(G * z - g);
        
        if Dynamic_res_norm <= PIPG_Tol && norm([delta_x;delta_u]) <= PIPG_Tol
            break;
        end
        p = p + 1;
    end

    disp(max(abs(z_opt-z)));
    if max(abs(z_opt-z)) >= 7e-4
        1;
    end
end
