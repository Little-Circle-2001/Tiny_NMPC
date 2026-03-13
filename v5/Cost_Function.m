function J = Cost_Function(state,input)
    StateDim = size(state,1);
    InputDim = size(input,1);
    N = size(input,2);

    Q = zeros(StateDim, StateDim, N+1);
    R = zeros(InputDim, InputDim, N);

    mu = 3.986e14;
    a = 6871393;
    n = sqrt(mu/a^3);
    
    x = state(1,:);
    y = state(2,:);
    z = state(3,:);
    vx = state(4,:);
    vy = state(5,:);
    vz = state(6,:);

    ux = input(1,:);
    uy = input(2,:);
    uz = input(3,:);
    
    for k = 1:N
        Q(:,:,k) = 1 * diag([1,1,1,1,1,1]);
        R(:,:,k) = 1 * eye(InputDim);
    end
    Q(:,:,N+1) = 10 * diag([1,1,1,1,1,1]);

    J = 0;
    for k = 1:N
        J = J + 1/2*state(:,k)'*Q(:,:,k)*state(:,k) + 1/2*input(:,k)'*R(:,:,k)*input(:,k);
    end
    J = J + 1/2*state(:,N+1)'*Q(:,:,N+1)*state(:,N+1);
end

