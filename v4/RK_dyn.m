function xk_next = RK_dyn(xk,uk,dt)
    %% step 1
    K1 = Quad_dyn(xk,uk);
    
    %% step 2
    K2 = Quad_dyn(xk+dt/2*K1,uk);

    %% step 3
    K3 = Quad_dyn(xk+dt/2*K2,uk);

    %% step 4
    K4 = Quad_dyn(xk+dt*K3,uk);

    xk_next = xk + dt/6*(K1+2*K2+2*K3+K4);
end