clc; clear;
close all;

monitor = trainingProgressMonitor;
monitor.Info = ["Iteration","LearningRate"];
monitor.Metrics = ["Loss","Residual"];
monitor.XLabel = 'Epoch';

%% 2. 生成训练数据（随机对称矩阵）
A = zeros(32,32,100);
Adata = zeros(32,32,100);
for i = 1:100
    X = randn(32); 
    A(:,:,i) = (X + X')/2;   % 对称化矩阵
    Adata(:,:,i) = A(:,:,i) / norm(A(:,:,i),'fro');
end

n = size(A,2);              
InputDim = n * n;            % 输入维度
hiddenSizes = [256 256];     % 两个隐藏层：256个和256个神经元
OutputDim = n + 1;           % 输出维度
SampleNum = size(A,3);       % 样本数


layers = [
featureInputLayer(InputDim, 'Name', 'input')

fullyConnectedLayer(hiddenSizes(1), 'Name', 'fc1')
reluLayer('Name', 'relu1')

fullyConnectedLayer(hiddenSizes(2), 'Name', 'fc2')
reluLayer('Name', 'relu2')

fullyConnectedLayer(OutputDim, 'Name', 'fc_out')    
];
net = dlnetwork(layers);

%% 训练参数
learnRate = 1e-3;
avgGrad = [];
avgSqGrad = [];
iteration = 0;
GradientDecayFactor = 0.9;
SquaredGradientDecayFactor = 0.999;
Epsilon = 1e-8;
Batch_size = 16;
alpha = 50.0;
kappa = 1.0;
alpha_min = 1e-4; 
kappa_min = 0;
alpha_decay = 0.8;
kappa_decay = 0.999;
numEpochs = 3000;

for epoch = 1 : numEpochs
    lossTotal = 0;
    resTotal = 0;
    for i = 1 : Batch_size: SampleNum
        iteration = iteration + 1;
        idx = i : min(i+Batch_size-1, SampleNum);

        % 求梯度并更新
        [loss, gradients, residualMean] = dlfeval(@NN_Loss, net, Adata(:,:,idx), alpha, kappa);
        lossTotal = lossTotal + extractdata(loss);
        resTotal = resTotal + extractdata(residualMean);
        [net.Learnables, avgGrad, avgSqGrad] = adamupdate(net.Learnables, ...
            gradients, avgGrad, avgSqGrad, iteration,...
            learnRate, GradientDecayFactor, SquaredGradientDecayFactor, Epsilon);
    end
    alpha = max(alpha_min,alpha*alpha_decay);
    kappa = max(kappa_min,kappa*kappa_decay);
    epoch_avg_loss = lossTotal / ceil(SampleNum/Batch_size);
    epoch_avg_res = resTotal / ceil(SampleNum/Batch_size);

    recordMetrics(monitor, epoch, ...
                  "Loss", loss, ...
                  "Residual", residualMean);
    updateInfo(monitor, ...
               "Iteration", iteration, ...
               "LearningRate", learnRate);


    if mod(epoch,100)==0
        fprintf('Epoch %d, Loss %.3e, Res %.3e, Alpha %.3e, Kappa %.3e\n', epoch, epoch_avg_loss, epoch_avg_res, alpha, kappa);
    end
end

save('TrainedNet.mat', 'net');
save('A.mat', "A");



