clear;clc;
syms alpha real;
x = [2;1];
delta_x = [1;1];
% f = @(alpha) (x + alpha * delta_x)' * (x + alpha * delta_x);
f = @(alpha) norm(x + alpha * delta_x);

y = f(alpha);
dy_sym = diff(y)
dy = matlabFunction(dy_sym)
subs(dy_sym,alpha,1)
dy(0)