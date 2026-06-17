function [alpha_star] = Line_Search_Method(f,px,pu,x0,u0,dt,c1,c2,MaxInteration)
    alpha = sym('alpha');
    
    phi_sym = f(x0+alpha*px, u0+alpha*pu, dt);
    der_phi_sym = diff(phi_sym,alpha);
%     phi_alpha = subs(phi_sym,{X,U},{x0,u0});
%     der_phi_alpha = subs(der_phi_sym,{X,U},{x0,u0});
%     phi_alpha = subs(subs(phi_sym, X, x0), U, u0);
%     der_phi_alpha = subs(subs(der_phi_sym, X, x0), U, u0);
    
    phi = matlabFunction(phi_sym);
    der_phi = matlabFunction(der_phi_sym);
    
    %% Line Search Algorithm
    alpha_star = Search_Alpha(phi,der_phi,c1,c2,MaxInteration);
end

