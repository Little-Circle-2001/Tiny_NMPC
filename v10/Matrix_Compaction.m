function [H, h, G, g] = Matrix_Compaction(Q, q, S, R, r, A, B, x0, gamma, N)
    H = [];    
    A_blk = [];
    B_blk = [];
    S_blk = [];
    
    n = size(A,1);
    m = size(B,2);
    G = zeros(N*n);
    g = zeros(N*n, 1);
    
    for k = 2 : N
        H = blkdiag(H, Q(:,:,k));
    end
    P = Q(:,:,N+1);
    H = blkdiag(H, P);
    for k = 1 : N 
        H = blkdiag(H, R(:,:,k));
    end

    h = zeros(size(H,2),1);
    for k = 1 : N
        h(n*(k-1)+1:n*k) = q(:,k+1);
        h(n*N+m*(k-1)+1:n*N+m*k) = r(:,k);
    end

    g(1:n) = A(:,:,1) * x0(:,1) + gamma(:,2);
    for k = 1 : N-1
        g(n*k+1:n*(k+1)) = gamma(:,k+2);
        A_blk = blkdiag(A_blk, -A(:,:,k+1));
        B_blk = blkdiag(B_blk, -B(:,:,k));
        S_blk = blkdiag(S_blk, S(:,:,k+1));
    end
    B_blk = blkdiag(B_blk, -B(:,:,k));
    G(n+1:end, 1:end-n) = A_blk;
    G = [G + eye(N*n), B_blk];

    H(1:(N-1)*n, N*n+m+1:end) = S_blk;
    H(N*n+m+1:end, 1:(N-1)*n) = S_blk';
end

