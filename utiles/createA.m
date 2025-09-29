function A = createA(X)
% 高效计算归一化余弦相似度矩阵
% 输入:
%   X - 数据矩阵 (fea x nSmp), 每列为一个样本
% 输出:
%   A - 归一化的相似度矩阵 (nSmp x nSmp)

[fea, nSmp] = size(X);

%% ========== 1. 计算余弦相似度矩阵 ==========
% 向量化计算余弦相似度 (替代原双重循环)
X_normalized = X ./ sqrt(sum(X.^2, 1)); % 按列归一化
AG = X_normalized' * X_normalized;      % 余弦相似度矩阵
AG = max(AG, 0);                       % 确保非负性(可选)
AG(1:nSmp+1:end) = 1;                 % 对角线置1

%% ========== 2. 对称归一化 ========== 
D = sum(AG, 2);                       % 节点度向量
D_inv_sqrt = spdiags(1./sqrt(D), 0, nSmp, nSmp); % 构造稀疏对角矩阵
A = D_inv_sqrt * AG * D_inv_sqrt;      % 归一化切割

%% 可选：转换为全矩阵（如果后续需要密集矩阵操作）
% A = full(A); 
end