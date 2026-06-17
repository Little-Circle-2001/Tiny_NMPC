function [Q, S, R, q, r, A, B, gamma] = Matrix_Component(H, h, G, g, StateDim, InputDim, N)
    Q = zeros(StateDim,StateDim,N+1);  
    S = zeros(StateDim,InputDim,N);
    R = zeros(InputDim,InputDim,N);
    q = zeros(StateDim,N+1);
    r = zeros(InputDim,N);
    A = zeros(StateDim,StateDim,N);
    B = zeros(StateDim,InputDim,N);
    gamma = zeros(StateDim,N+1);
    for k = 2 : N
        Q(:,:,k) = H(StateDim*(k-2)+1 : StateDim*(k-1),...
                     StateDim*(k-2)+1 : StateDim*(k-1));
        S(:,:,k) = H(StateDim*(k-2)+1 : StateDim*(k-1),...
                     StateDim*N+InputDim*(k-1)+1 : StateDim*N+InputDim*k);
        R(:,:,k-1) = H(StateDim*N+InputDim*(k-2)+1 : StateDim*N+InputDim*(k-1),...
                       StateDim*N+InputDim*(k-2)+1 : StateDim*N+InputDim*(k-1));

        q(:,k) = h(StateDim*(k-2)+1 : StateDim*(k-1));
        r(:,k-1) = h(StateDim*N+InputDim*(k-2)+1 : StateDim*N+InputDim*(k-1));

        A(:,:,k) = -1*G(StateDim*(k-1)+1 : StateDim*k,...
                        StateDim*(k-2)+1 : StateDim*(k-1));
        B(:,:,k-1) = -1 * G(StateDim*(k-2)+1 : StateDim*(k-1),...
                            StateDim*N+InputDim*(k-2)+1 : StateDim*N+InputDim*(k-1));
        gamma(:,k) = g(StateDim*(k-2)+1 : StateDim*(k-1));
    end
    Q(:,:,N+1) = H(StateDim*(N-1)+1 : StateDim*N,...
                   StateDim*(N-1)+1 : StateDim*N);
    R(:,:,N) = H(StateDim*N+InputDim*(N-1)+1 : StateDim*N+InputDim*N,...
                   StateDim*N+InputDim*(N-1)+1 : StateDim*N+InputDim*N);

    q(:,N+1) = h(StateDim*(N-1)+1 : StateDim*N);
    r(:,N) = h(StateDim*N+InputDim*(N-1)+1 : StateDim*N+InputDim*N);

    B(:,:,N) = -1 * G(StateDim*(N-1)+1 : StateDim*N,...
                      StateDim*N+InputDim*(N-1)+1 : StateDim*N+InputDim*N);
    gamma(:,N+1) = g(StateDim*(N-1)+1 : StateDim*N);
end

