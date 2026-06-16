function [x_iter, u_iter, J_iter, SCP_Iter_Num] = PIPG_Based_SCP(x0, u0, target, ...
                                                                 Obs_Pos, Obs_r, ...
                                                                 SCP_MaxIter, SCP_Tol, ...
                                                                 eta, rho, selfsigma, ADMM_MaxIter, ADMM_abs_eps, ADMM_rel_eps, ...
                                                                 rho_update_inteval, adaptive_rho_tolerance, ...
                                                                 PIPG_MaxIter, PIPG_Tol, omega, ...
                                                                 lambda0, v0, mu0)
    StateDim = size(x0,1);
    InputDim = size(u0,1);
    N = size(u0,2);
    dt = 6;
    
    delta_x0 = zeros(StateDim, N+1);
    delta_u0 = zeros(InputDim, N);
    if nargin == 18
        [~, ~, ~, ~, ~, ~, ~, ~, D, ~, ~, ~, ~] = Make_Model(x0, u0, target, Obs_Pos, Obs_r);
        lambda0 = zeros(StateDim,N+1);
        mu0 = cell(1,N+1);
        v0 = cell(1,N+1);
        for i = 1:N+1
            mu0{i} = zeros(size(D{i},1), 1);
            v0{i} = zeros(size(D{i},1), 1);
        end
    end
    %% Intialize
    % Corresponding to Line 2-4 of Algorithm 2 in the manuscripts
    SCP_Iter_Num = 0;
    J0 = Merit_Function(x0,u0,dt);
    x1 = x0;
    u1 = u0;
    J1 = J0;
    x_iter = x1;
    u_iter = u1;
    J_iter = J1;
    
    tic;
    %% Start SCP
    delta_J = 10;
    delta_x = 10;
    delta_u = 10;
    while (SCP_Iter_Num < SCP_MaxIter) && ((abs(delta_J) > SCP_Tol) && abs(delta_x) > SCP_Tol && abs(delta_u) > SCP_Tol )
        % Calculate H,h,G,g via Eq.(15)-(16)
        % Corresponding to Line 6 of Algorithm 2 in the manuscripts
        [Q, q, S, R, r, A, B, gamma, D, E, c, PIPG_z_lower, PIPG_z_upper] = Make_Model(x0, u0, target, Obs_Pos, Obs_r);
%         for k = 1:N
%             Q(:,:,k) = Q(:,:,k) + rho*D{k}'*D{k} + selfsigma*eye(StateDim);       % Update \tilde{Q_k}
%             R(:,:,k) = R(:,:,k) + rho*E{k}'*E{k} + selfsigma*eye(InputDim);       % Update \tilde{R_k}
%             S(:,:,k) = S(:,:,k) + rho*D{k}'*E{k};                                 % Update \tilde{S_k}
%         end        
%         Q(:,:,N+1) = Q(:,:,N+1) + rho*D{N+1}'*D{N+1} + selfsigma*eye(StateDim);   % Update \tilde{Q_N}
%         [H, h, G, g] = Matrix_Compaction(Q, q, S, R, r, A, B, delta_x0, gamma, N);
% 
%         % Calculate α and β
%         % Corresponding to Line 7 of Algorithm 2
%         nu = max(eig(H));
%         sigma = max(eig(G'*G));
%         alpha = 2/(nu+sqrt(nu^2+4*omega*sigma));
%         beta = omega * alpha;

        % ADMM求解
        [delta_x0, delta_u0, lambda1, mu1] = ADMM_Solver(Q, q, S, R, r, A, B, gamma, D, E, c, ...                                     
                                                         delta_x0, delta_u0, lambda0, v0, mu0, ...
                                                         eta, rho, selfsigma, ADMM_MaxIter, ADMM_abs_eps, ADMM_rel_eps, ...
                                                         rho_update_inteval, adaptive_rho_tolerance, ...
                                                         PIPG_MaxIter, PIPG_Tol, omega, PIPG_z_lower, PIPG_z_upper);
        alpha_star = 1;
        Use_Linear_Search = 1;
        if Use_Linear_Search == 1
            MaxInteration = 10;
            for i = 1:MaxInteration
                if i == MaxInteration
%                     print("No improvement")
                    alpha_star = 0.5*alpha_star;
                    break;
                end
                x_try = x0+alpha_star * delta_x0;
                u_try = u0+alpha_star * delta_u0;
                if Merit_Function(x_try,u_try,dt) > Merit_Function(x0,u0,dt)
                    alpha_star = 0.5 * alpha_star;
                else
                    break;
                end
            end
        else
            alpha_star = 1;
        end

        % State and Input variables Update
        % Corresponding to Line 15-16 of Algorithm 2
        x1 = x0 + alpha_star * delta_x0;
        u1 = u0 + alpha_star * delta_u0;
        J1 = Merit_Function(x1,u1,dt);
        x_iter = cat(3, x_iter, x1);
        u_iter = cat(3, u_iter, u1);
        J_iter = cat(2, J_iter, J1);
        delta_J = J0 - J1;
        delta_x = norm(x0 - x1);
        delta_u = norm(u0 - u1);

        % Update intial variables
        x0 = x1;
        u0 = u1;
        J0 = J1;
        lambda0 = lambda1;
        mu0 = mu1;

        % j = j + 1
        SCP_Iter_Num = SCP_Iter_Num + 1;
    end

end