%% zoom函数
%
%%
function [alpha_star] = zoom(phi,der_phi,c1,c2,iter_num,alpha_lo,alpha_hi)
    for i = 1:iter_num
%         alpha_j = -(der_phi(alpha_lo)*alpha_hi^2)/(2*(phi(alpha_hi)-phi(alpha_lo)-der_phi(alpha_lo)*alpha_hi));
        alpha_j = (alpha_hi + alpha_lo)/2;
        value = phi(alpha_j);
        if (value>phi(0)+c1*alpha_j*der_phi(0)) || value >= phi(alpha_lo)
            alpha_hi = alpha_j;
        else
            gradient_value = der_phi(alpha_j);
            if abs(gradient_value)<=-c2*der_phi(0)
                alpha_star = alpha_j;
                return;
            end
            if gradient_value*(alpha_hi-alpha_lo)>=0
                alpha_hi = alpha_lo;
            end
            alpha_lo = alpha_j;
        end
        alpha_star = alpha_j;
    end
end

