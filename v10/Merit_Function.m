function J = Merit_Function(state,input,dt)
    StateDim = size(state,1);
    InputDim = size(input,1);
    N = size(input,2);

    Q = zeros(StateDim, StateDim, N+1);
    R = zeros(InputDim, InputDim, N);

    mu = 3.986e5;
    a = 6871.393;
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
        Q(:,:,k) = 1 * diag([1e-5,1e-5,1e-5,1,1,1]);
        R(:,:,k) = eye(InputDim);
    end
    Q(:,:,N+1) = diag([1,1,1,1,1,1]);

    J = 0;
    for k = 1:N
        J = J + 1/2*state(:,k)'*Q(:,:,k)*state(:,k) + 1/2*input(:,k)'*R(:,:,k)*input(:,k);
    end
    J = J + 1/2*state(:,N+1)'*Q(:,:,N+1)*state(:,N+1);

%     A = [    0      0      0      1      0      0; ...
%              0      0      0      0      1      0; ...
%              0      0      0      0      0      1; ...
%          3*n^2      0      0      0    2*n      0; ...
%              0      0      0   -2*n      0      0; ...
%              0      0   -n^2      0      0      0; ...
%         ] * dt;
%     B = [0      0      0; ...
%          0      0      0; ...
%          0      0      0; ...
%          1      0      0; ...
%          0      1      0; ...
%          0      0      1
%          ] * dt;

    for k = 1:N
        J = J + 1000 * norm(state(:,k+1) - state(:,k) - [vx(k); ...
                                                         vy(k); ...
                                                         vz(k); ...
                                3*n^2*x(k) + 2*n*vy(k) + ux(k); ...
                                            -2*n*vx(k) + uy(k); ...
                                             -n^2*z(k) + uz(k)] * dt);
    end

%     for 
end

