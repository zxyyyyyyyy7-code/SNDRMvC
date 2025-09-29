clc;
clear;
rng('default');
rng(5489);
% rng(2021);
addpath('datasets/');
options = [];
% 数据集配置列表
datasets = {
   struct('name', 'MSRC', 'file', 'MSRC.mat'),
    % 继续添加更多数据集...
};

 
% 主循环：每个数据集依次处理
if ~exist('result', 'dir')
    mkdir('result');
end
for d = 1:length(datasets)
    ds = datasets{d};
     diary(fullfile('result', ['output_' ds.name '.txt']));
    fprintf('>>> Running on dataset: %s\n', ds.name);

    % =============== 每个数据集单独处理 ===============
    switch ds.name  
        case 'MSRC'
            load(ds.file); % fea, gt
            data = cellfun(@transpose, X, 'UniformOutput', false);
            gnd = Y;
            options.eta=[0.001];
            options.etaH=[0.001];
            options.beta=[0.001];
            options.K=[10];   %低维维度
  
        otherwise
            error(['未知数据集: ' ds.name]);
    end
   

    %% ————————————————————options初始化————————————————————————————————————————
        n=size(gnd);
       
        options.maxIter = 100;
        options.error = 1e-6;
        options.nRepeat = 10;
        options.minIter = 30;
        options.meanFitRatio = 0.1;
        options.rounds = 10;
        options.Gaplpha=1; %Graph regularisation parameter
        options.WeightMode='Binary';
        %————————计算相似性矩阵参数————————
        options.graph_type='sparse';%稀疏图
        %options.graph_type='sparse';%全连接图
        options.similarity_type='gaussian';%高斯相似性
        %options.similarity_type='inner_product';%内积相似性
        options.kk = floor(log2(n)) + 1;
        options.nn = 7;
        %——————————————————————————————
        options.nClass=length(unique(gnd));
        options.kmeans = 1;

    tic;
    % 归一化 + 构图
    A = cell(1, length(data));
    for i = 1:length(data)
        data{i} = NormalizeFea1(data{i});
        A{i} = createA1(data{i}, options.graph_type, options.similarity_type, options.kk, options.nn);
    end

    % 算法运行
    for eta = options.eta
        for beta = options.beta
            for K = options.K
                for etaH = options.etaH
                    fprintf('Eta: %.6f, EtaH: %.6f, Beta: %.6f, K: %.2f\n', eta, etaH, beta, K);
                     [LABEL] = SNDRMvC(A, K , options,gnd,eta,etaH,beta);
                     time = toc;
                    res=Clustering8Measure(LABEL{30},gnd);
                    disp(res)  
                    fprintf('Time for : %.4f seconds\n', time);
                end
            end
        end
    end

    diary off;
end
