function [LABEL] = SNDRMvC(X, K , options, gnd,eta,etaH,beta)
rng('default')
rng(5489)

% 初始化参数
C = options.nClass; % 类别数目
% 设置迭代轮数
Round = 0;
maxRound=30;
% 获取视图数量
nView = numel(X);
lumdav = 1/nView;

obj_NMF1 = cell(1, nView);
obj_NMF2 = cell(1, nView);
obj_NMF3 = cell(1, nView);
obj_NMF4 = cell(1, nView);
obj_NMF = cell(1, nView);
Ux = cell(1, nView);
dX=cell(1, nView);


% 初始化 H 为单元格数组，每个单元格存储一个视图的低维表示
H = cell(1, nView);
% 初始化 m 为单元格数组，每个单元格存储一个视图的偏置项
m = cell(1, nView);
% 初始化 h 和 m
for viewIdx = 1:nView
    [mFea, nSmp] = size(X{viewIdx});
    H{viewIdx} = logsig(rand(mFea, K)); % 初始化 h
    m{viewIdx} = logsig(rand(mFea, 1)); % 初始化 m 
end
% centroidH = logsig(rand(mFea, C));
centroidH=H{1};


% 初始化 ACC 和 NMI 结果数组
acc = zeros(1, 30);
NMI = zeros(1, 30);
LABEL = cell(1, 30); % 存储每轮的标签

% 迭代更新 h 和 m
while Round < maxRound
    Round = Round + 1;
   
    % 循环遍历每个视图
    for viewIdx = 1:nView

        %取出当前视图
        current_X=X{viewIdx};
        current_H=H{viewIdx};
        current_m=m{viewIdx};

        % 获取当前视图的尺寸
        [mFea, nSmp] = size(current_X);

        % 非线性约束 SNMF (NNMF)
        for i = 1:mFea
            for j = 1:nSmp
                sum1 = 0;

                for kk = 1:K
                    sum1 = sum1 + 1/(1 + exp(-current_H(i, kk))) * 1/(1 + exp(-current_H(j, kk)));
                end
                temp = current_X(i, j) - sum1 - 1/(1 + exp(-current_m(i)));

                % 更新偏置项 m
                C1 = 1/(1 + exp(-current_m(i)));
                current_m(i) = current_m(i) - eta * (temp * (-1) * C1 * (1 - C1));

                % 更新低维矩阵 h 的元素
                for kk = 1:K
                    A = 1/(1 + exp(-current_H(i, kk)));
                    CH = 1/(1 + exp(-centroidH(i, kk)));
                    current_H(i, kk) = current_H(i, kk) - eta * (temp * (-1)*1/(1+exp(-current_H(j, kk))) + lumdav * (A - CH) + beta * A )* (1 - A) * A;
                end
            end
        end

       
         % 计算低维矩阵 h 的值
        H{viewIdx} = current_H;
        % 计算偏置 m 的值
        m{viewIdx} = current_m;
        
    end
 % 计算共识矩阵的值
 for viewIdx = 1:nView
        %取出当前视图21
        current_H=H{viewIdx}; 
        for i = 1:mFea
            for kk = 1:K
                A = 1/(1 + exp(-current_H(i, kk)));
                CA = 1/(1 + exp(-centroidH(i, kk)));
                centroidH(i, kk) = centroidH(i, kk) - etaH * lumdav * ((-1) *(A-CA))*CA*(1-CA);
            end
        end
 end

 for viewIdx = 1:nView
    H{viewIdx} = logsig(H{viewIdx});
    m{viewIdx} = logsig(m{viewIdx});
 end
 centroidH = logsig(centroidH);

    obj_NMFall(Round) = 0;
    O1=0;
    O2=0;
    O3=0;
 for viewIdx = 1:nView
        current_H=H{viewIdx};
        current_m=m{viewIdx};
        Ux{viewIdx} = [current_H, current_m];
        dX{viewIdx} = Ux{viewIdx}*(Ux{viewIdx})'-X{viewIdx};
        obj_NMF1{viewIdx}(Round) = sum(sum(dX{viewIdx}.^2));
        obj_NMF2{viewIdx}(Round) = beta*sum(sum(current_H.^2));
        obj_NMF3{viewIdx}(Round) = lumdav*sum(sum(current_H-centroidH).^2);
        obj_NMF4{viewIdx}(Round)=obj_NMF1{viewIdx}(Round)+obj_NMF2{viewIdx}(Round)+obj_NMF3{viewIdx}(Round);
        obj_NMF{viewIdx}(Round)=sqrt(obj_NMF4{viewIdx}(Round));
        obj_NMFall(Round) = obj_NMFall(Round)+obj_NMF{viewIdx}(Round);
        O1=O1+obj_NMF1{viewIdx}(Round);
        O2=O2+obj_NMF2{viewIdx}(Round);
        O3=O3+obj_NMF3{viewIdx}(Round);  
 end
  disp(['obj_NMF1：',  num2str(O1)]);
    disp(['obj_NMF2：',   num2str(O2)]);
    disp(['obj_NMF3: ',   num2str(O3)]);
 
    % 获得聚类标签
    rng(5489)
    label = litekmeans(centroidH, C, 'Replicates', 20);

 % 确保标签向量是列向量
    gnd   = gnd(:);
    label = label(:);

    % 计算最佳映射（用于 ACC）
    idx22 = bestMap(gnd, label);

    % 计算 ACC
    LABEL{Round} = idx22;
    acc(Round) = sum(gnd == idx22) / length(gnd);

    % ---- 计算 NMI 前做标签映射到 1..K ----
    [~, ~, gnd_tmp]   = unique(gnd,   'stable');   % 与 gnd 等长
    [~, ~, label_tmp] = unique(label, 'stable');   % 与 label 等长

    % 检查异常
    if length(gnd_tmp) ~= length(label_tmp)
        error('标签长度不一致！gnd=%d, label=%d', length(gnd_tmp), length(label_tmp));
    end
    if any(isnan(gnd_tmp)) || any(isnan(label_tmp))
        error('标签中含 NaN，请检查数据。');
    end

 % 计算 NMI
    NMI(Round) = MutualInfo(gnd_tmp, label_tmp);
    disp(['Number of Rounds:', num2str(Round)]);
    disp(['Object Function Value: ',num2str(obj_NMFall(Round))]);
    disp(['Acc:', num2str(acc(Round))]);
    disp(['NMI:', num2str(NMI(Round))]);
    disp('-------------------------------------------------------');
end
end


