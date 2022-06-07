function [simulacao]  = simulate_markov_chain_finalDefined(P, n, initial_facies, final_facies, nsims, prob_map)
                        
if nargin<4
    nsims=1;
end

n_facies = size(P,2);

P_T = transpose_trasition_matrix(P);

simulacao = zeros(n,nsims);

simulacao(1,:) = initial_facies;
simulacao(end,:) = final_facies;

if nargin<5
    for j=1:nsims
        for i = 2:n-1   
        facies = find( rand < cumsum(P(simulacao(i-1,j),:)));
        simulacao(i,j) = facies(1);    
        end
    end
else
   for j=1:nsims
       for i = 2:n-2
       probabilities = P(simulacao(i-1,j),:).*reshape(prob_map(i,1,:),1,n_facies);
       probabilities = probabilities./sum(probabilities);
       facies = find( rand < cumsum( probabilities ));
       simulacao(i,j) = facies(1);    
       end
       i=n-1;
       probabilities = P(simulacao(i-1,j),:).*P_T(simulacao(i+1,j),:).*reshape(prob_map(i,1,:),1,n_facies);
       probabilities = probabilities./sum(probabilities);
       facies = find( rand < cumsum( probabilities ));
       simulacao(i,j) = facies(1);    
   end
end
