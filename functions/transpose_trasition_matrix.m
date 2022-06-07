function [P_upward] = transpose_trasition_matrix(P)

n_facies = size(P,1);

prop = P^100;
prop = prop(1,:);

P_upward = P';
P_upward = P_upward.*repmat(prop,n_facies,1);
P_upward = P_upward./repmat(sum(P_upward,2),1,n_facies);

if sum(sum(isnan(P_upward))) >0
    P_upward = ones(n_facies,n_facies)/n_facies;
    %P_upward = P_upward./repmat(sum(P_upward,2),1,n_facies);
end