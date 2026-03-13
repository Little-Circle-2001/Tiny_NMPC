%% 寻找满足Wolfe Condition的步长的线搜索方法
% 
%% Line_Search_Algorithm 
function [alpha_star] = Search_Alpha(phi,der_phi,c1,c2,iter_num)
    % 定义初值
    alpha0 = 0;
    alpha_max = 2;
    alpha = 0.3;
    value0 = phi(0);
    for i = 1:iter_num       
        value = phi(alpha);
        if value > phi(0)+c1*alpha*der_phi(0) || (value>=value0 && i>1)
            alpha_star = zoom(phi,der_phi,c1,c2,iter_num,alpha0,alpha);
            break;
        end
        gradient_value = der_phi(alpha);
        if abs(gradient_value)<=abs(c2*der_phi(0))
            alpha_star = alpha;
            break;
        end
        if gradient_value>=0
            alpha_star = zoom(phi,der_phi,c1,c2,iter_num,alpha,alpha0);
            break;
        end
        alpha0 = alpha;
        value0 = value;
        alpha = min(1.5*alpha,alpha_max);
    end

end

