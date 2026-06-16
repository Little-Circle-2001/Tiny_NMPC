function J = Cost_Function(x,u)
    StateDim = size(x,1);
    InputDim = size(u,1);
    N = size(u,2);

    Q = zeros(StateDim, StateDim, N+1);
    R = zeros(InputDim, InputDim, N);
    
    for k = 1:N
        Q(:,:,k) = 10 * diag([2,2,2,0.5,0.5,0.5,0.1,0.1,0.1,0.1,0.1,0.1]);
        R(:,:,k) = 5 * eye(InputDim);
    end
    Q(:,:,N+1) = 20 * diag([2,2,2,0.5,0.5,0.5,0.1,0.1,0.1,0.1,0.1,0.1]);

    J = 0;
    for k = 1:N
            J = J + 1/2*x(:,k)'*Q(:,:,k)*x(:,k) + 1/2*u(:,k)'*R(:,:,k)*u(:,k);
    end
    J = J + 1/2*x(:,N+1)'*Q(:,:,N+1)*x(:,N+1);
end

