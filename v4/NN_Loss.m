function [loss, gradients, residualMean] = NN_Loss(net, Adata, alpha, kappa)
    num = size(Adata,3);
    batch_loss = dlarray(0.0);
    batch_res = dlarray(0.0);
    for i = 1:num
        A = Adata(:,:,i);
        Aflat = A(:);
        dlA = dlarray(Aflat,'CB');
        dlout = forward(net,dlA);
    
        lambda = dlout(1);
        v = dlout(2:end);
    
        % 残差项
        res = A*v - lambda*v;
        loss_res = sum(res.^2);
        
        % 归一化项
        loss_norm = (sum(v.^2) - 1)^2;
    
        % 线性项
        loss_sqa = 0.5 * v'*A*v;
    
        % 总损失
        s = -1; % -1 -> 最大特征值方向
        batch_loss = batch_loss + s * alpha * loss_sqa + loss_res + kappa * loss_norm;
        batch_res = batch_res + loss_res;
    end
    loss = batch_loss / num;
    residualMean = batch_res / num;
    gradients = dlgradient(loss, net.Learnables);
end