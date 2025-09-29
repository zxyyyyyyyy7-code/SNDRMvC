function [A] = createA1(X,graph_type,similarity_type,kk,nn)
D = dist2(X, X);
if strcmp(graph_type, 'full') & strcmp(similarity_type, 'gaussian')
    A = scale_dist3(D, nn);
elseif strcmp(graph_type, 'full') & strcmp(similarity_type, 'inner_product')
    A = X * X';
elseif strcmp(graph_type, 'sparse') & strcmp(similarity_type, 'gaussian')
    A = scale_dist3_knn(D, nn, kk, true);
else % graph_type == 'sparse' & similarity_type == 'inner_product'
    Xnorm = X';
    d = 1./sqrt(sum(Xnorm.^2));
    Xnorm = bsxfun(@times, Xnorm, d);
    A = inner_product_knn(D, Xnorm, knn, true);
    clear Xnorm, d;
end
clear D;


dd = 1 ./ sum(A);
dd = sqrt(dd);
A = bsxfun(@times, A, dd);
A = A';
A = bsxfun(@times, A, dd);
clear dd;

A = (A + A') / 2;
end

