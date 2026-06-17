function xk_next = RK_dyn(xk,uk,dt,n)
    %% step 1
    K1 = CW_dyn(xk,uk,n);
    
    %% step 2
    K2 = CW_dyn(xk+dt/2*K1,uk,n);

    %% step 3
    K3 = CW_dyn(xk+dt/2*K2,uk,n);

    %% step 4
    K4 = CW_dyn(xk+dt*K3,uk,n);

    xk_next = xk + dt/6*(K1+2*K2+2*K3+K4);
end