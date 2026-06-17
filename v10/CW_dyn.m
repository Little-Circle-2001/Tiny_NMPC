function dx = CW_dyn(xk,uk,n)
    dx =zeros(6,1);
    x = xk(1);
    y = xk(2);
    z = xk(3);
    vx = xk(4);
    vy = xk(5);
    vz = xk(6);
    ux = uk(1);
    uy = uk(2);
    uz = uk(3);

    dx(1) = vx;
    dx(2) = vy; 
    dx(3) = vz;
    dx(4) = 3*n^2*x + 2*n*vy + ux;
    dx(5) = -2*n*vx + uy;
    dx(6) = -n^2*z + uz;                       
end