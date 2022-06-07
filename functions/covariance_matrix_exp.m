function [covar] = covariance_matrix_exp(sgm2,L,order)

sgm = sqrt(sgm2);
I = length(sgm);

[X,Y] = meshgrid([1:I],[1:I]);

covar = exp(-(abs(X(:)-Y(:))/L).^order);
covar = reshape(covar,I,I);

SGM = diag(sqrt(sgm2));

covar = SGM *covar * SGM ;
