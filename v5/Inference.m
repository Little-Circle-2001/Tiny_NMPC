clear; clc;
close all;
load('A.mat', 'A');
load("TrainedNet.mat", "net")
i = 100;
Apred = A(:,:,i);
Apredflat = dlarray(Apred(:), "CB");
Ypred = predict(net, Apredflat);
lambda = Ypred(1);
v = Ypred(2:end);
res = Apred*v - lambda*v;
disp(v'*v);
disp(res);


[V, D] = eig(Apred);
lambdas_Act = diag(D);               % 提取特征值为列向量
[lambda_Act, idx] = max(lambdas_Act); % 找到最大特征值
v_Act = V(:, idx);               % 对应特征向量
res_Act = Apred*v_Act - lambda_Act*v_Act;
disp(res_Act);

