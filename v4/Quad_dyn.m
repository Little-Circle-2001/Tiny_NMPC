function dx = Quad_dyn(xk,uk)
    dx =zeros(12,1); m = 2.69; g = 9.81;
    Ixx = 0.015;Iyy = 0.015;Izz = 0.0245;

    pN = xk(1);pE = xk(2);pD = xk(3);vN = xk(4);vE = xk(5);vD = xk(6);
    phi = xk(7);theta = xk(8);psi = xk(9);omega_x = xk(10);omega_y = xk(11);omega_z = xk(12);
    F = uk(1);taux = uk(2);tauy = uk(3);tauz = uk(4);

    dx(1) = vN;
    dx(2) = vE;
    dx(3) = vD;
    dx(4) = -(sin(theta)*cos(phi)*cos(psi)+sin(psi)*sin(phi))*F/m;
    dx(5) = -(sin(theta)*sin(psi)*cos(phi)-sin(phi)*cos(psi))*F/m;
    dx(6) = (m*g-F*(cos(theta)*cos(phi)))/m;
    dx(7) = omega_x+omega_y*sin(phi)*tan(theta)+omega_z*cos(phi)*tan(theta);
    dx(8) = omega_y*cos(phi)-omega_z*sin(phi);
    dx(9) = 1/cos(theta)*(omega_y*sin(phi)+omega_z*cos(phi));
    dx(10) = ((Iyy-Izz)*omega_y*omega_z+taux)/Ixx;
    dx(11) = ((Izz-Ixx)*omega_x*omega_z+tauy)/Iyy;
    dx(12) = ((Ixx-Iyy)*omega_x*omega_y+tauz)/Izz;   
end