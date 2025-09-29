function fea = NormalizeFea1(fea, row, pnorm)
if nargin < 2, row = 1; end
if nargin < 3, pnorm = 2; end

if row
    for i = 1:size(fea,1)
        norm_i = max(1e-14, sum(abs(fea(i,:)).^pnorm)^(1/pnorm));
        fea(i,:) = fea(i,:) / norm_i;
    end
else
    for i = 1:size(fea,2)
        norm_i = max(1e-14, sum(abs(fea(:,i)).^pnorm)^(1/pnorm));
        fea(:,i) = fea(:,i) / norm_i;
    end
end
end
