function [simulation]  = simulate_markov_chain(P, n, initial_facies, nsims, prob_map)
% output:   simulation - 1D facies simulations (n X nsims)-matrix, each column is a 1D simulation
% input:    P - Transition matrix
%           n - Size of the simulation
%           initial_facies - Initial facies of the chain to start the simulation
%           nsims - Number of simulations
%           prob_map (optional) - Pointwise prior probability, usually comes from the Bayesian inference/classification

if nargin<4
    nsims=1;
end

n_facies = size(P,2);

simulation = zeros(n,nsims);

simulation(1,:) = initial_facies;

if nargin<5
    for j=1:nsims
        for i = 2:n   
        facies = find( rand < cumsum(P(simulation(i-1,j),:)));
        simulation(i,j) = facies(1);    
        end
    end
else
   for j=1:nsims
       for i = 2:n   
       probabilities = P(simulation(i-1,j),:).*reshape(prob_map(i,1,:),1,n_facies );
       probabilities = probabilities./sum(probabilities);
       facies = find( rand < cumsum( probabilities ));
       simulation(i,j) = facies(1);    
       end
   end
end