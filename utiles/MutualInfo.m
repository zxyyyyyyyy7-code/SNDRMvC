function MIhat = MutualInfo(L1, L2)
%   mutual information
%
%   version 2.0 --May/2007
%   version 1.0 --November/2003
%
%   Written by Deng Cai (dengcai AT gmail.com)
%===========    
L1 = L1(:);
L2 = L2(:);

% 检查输入数据是否包含无效值
if any(isnan(L1)) || any(isnan(L2))
    error('L1 or L2 contains NaN values');
end
if any(isinf(L1)) || any(isinf(L2))
    error('L1 or L2 contains Inf values');
end

if size(L1) ~= size(L2)
    error('size(L1) must == size(L2)');
end

Label = unique(L1);
nClass = length(Label);

Label2 = unique(L2);
nClass2 = length(Label2);

% 动态平滑
if nClass2 < nClass
    missingLabels = setdiff(Label, unique(L2));
    L1 = [L1; missingLabels];
    L2 = [L2; missingLabels];
elseif nClass2 > nClass
    missingLabels = setdiff(Label2, unique(L1));
    L1 = [L1; missingLabels];
    L2 = [L2; missingLabels];
end

G = zeros(nClass);
for i = 1:nClass
    for j = 1:nClass
        G(i, j) = sum(L1 == Label(i) & L2 == Label(j));
    end
end
sumG = sum(G(:));

% 添加小常数进行平滑
epsilon = 1e-10;
P1 = sum(G, 2); P1 = (P1 + epsilon) / (sumG + nClass * epsilon);
P2 = sum(G, 1); P2 = (P2 + epsilon) / (sumG + nClass * epsilon);

if sum(P1 == 0) > 0 || sum(P2 == 0) > 0
    % 如果仍然存在零值，抛出错误
    error('Smooth fail!');
else
    H1 = sum(-P1 .* log2(P1));
    H2 = sum(-P2 .* log2(P2));
    P12 = G / sumG;
    PPP = P12 ./ repmat(P2, nClass, 1) ./ repmat(P1, 1, nClass);
    PPP(abs(PPP) < 1e-12) = 1;
    MI = sum(P12(:) .* log2(PPP(:)));
    MIhat = MI / max(H1, H2);
    %%%%%%%%%%%%%   why complex ?       %%%%%%%%
    MIhat = real(MIhat);
end