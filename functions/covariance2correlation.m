function [nu] = covariance2correlation(C)

I = size(C,1);

stds = sqrt(diag(C));

nu = C./repmat(stds',I,1);
nu = nu./repmat(stds,1,I);

