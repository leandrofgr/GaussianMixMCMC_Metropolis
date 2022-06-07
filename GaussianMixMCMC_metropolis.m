function [ INVERSION ] = GaussianMixMCMC_metropolis(real_seismic, theta, signal2noise, wavelet, PRIOR, n_it, prob_map, P, facies_initial)

tic

%% burn in period
it_conv = min([10000 round(n_it*0.1)]);

n_angles = length(theta);
I = size(real_seismic,1);
n_facies = length(PRIOR);

real_seismic  = real_seismic(1:end-1,:);

G = elasticForwardModel(exp(PRIOR(1).MU(1))*ones(I,1), exp(PRIOR(1).MU(2))*ones(I,1), wavelet, theta);

%% Seismic noise covariance matrix 
d = reshape(real_seismic,size(real_seismic,2)*size(real_seismic,1),1);
var_sismica_final = var(real_seismic)./(signal2noise);
var_d = repmat(var_sismica_final,size(d,1)/n_angles,1);
C_d = diag(var_d(:));

I_d = size(d,1);

%% INITIAL CONFIGURATION - RANDOM MODEL
if nargin < 9
    facies_initial = simulate_markov_chain(P, I, randi(n_facies), 1, prob_map);
end

%% Compute the prior distribution of elastic properties along the traces given de facies configuration
facies_samples(:,1) = facies_initial;
if length(fieldnames(PRIOR))==2
    [mu_mlpi_previous , C_mlpi_previous ] = construct_prior_facies(PRIOR,facies_initial,2.5);
else
    for facie=1:length(PRIOR)
        trendVp(:,facie) = PRIOR(facie).VP.trend;
        trendVs(:,facie) = PRIOR(facie).VS.trend;
        trendRhob(:,facie) = PRIOR(facie).RHOB.trend;
    end
    [mu_mlpi_previous , C_mlpi_previous ] = construct_prior_trend_correlated_v2(PRIOR, trendVp, trendVs, trendRhob, facies_initial, 2.5);
end

%% MCMC METHOD
model_samples = zeros(3*I,n_it);
sy_seisimic = zeros(I_d,n_it);
C_samples = zeros(3*I,n_it);
log_likelihood = zeros(1,n_it);
log_likelihood(1) = nan;
OPERATOR = C_mlpi_previous*G'*inv(G*C_mlpi_previous*G' + C_d);
model_samples(:,1) = mu_mlpi_previous + OPERATOR*(d - G*mu_mlpi_previous);
C_samples (:,1) = diag(C_mlpi_previous - OPERATOR*G*C_mlpi_previous);
for it = 2:n_it
    100*it/n_it
    
    %% From the previous facies configuration, drawn a facies samples from the prior
    window_size = round(4 + rand*6);
    window_size = window_size - 1;
    
    index_sort = 1 + round(rand*(I-window_size-1)) + window_size/2;
    indexes_sort = index_sort-window_size/2:window_size/2+index_sort;
    indexes_sort_not = [1:(indexes_sort(1)-1) (indexes_sort(end)+1):I ];
    
    if indexes_sort(1) == 1 || indexes_sort(end) == I
        if indexes_sort(1) == 1
            initial_facies = randi(n_facies);
            final_facies = facies_samples(indexes_sort(end)+1,it-1);
        else
            initial_facies = facies_samples(indexes_sort(1)-1,it-1);
            final_facies = randi(n_facies);
        end
    else
        initial_facies = facies_samples(indexes_sort(1),it-1);
        final_facies = facies_samples(indexes_sort(end),it-1);
    end
    
    prob_map_window = prob_map(indexes_sort,:,:);
    
    facies_proposal_window = simulate_markov_chain_finalDefined(P, window_size+1, initial_facies, final_facies, 1, prob_map_window);
    facies_proposal = facies_samples(:,it-1);
    facies_proposal(indexes_sort) = facies_proposal_window ;
    
    %% Construct prior mean an covariance of the the prior p(m|pi)
    if length(fieldnames(PRIOR))==2
        [mu_mlpi, C_mlpi] = construct_prior_facies(PRIOR,facies_proposal,2.5);
    else
        [mu_mlpi, C_mlpi] = construct_prior_trend_correlated_v2(PRIOR, trendVp, trendVs, trendRhob, facies_proposal,2.5);
    end
    
    %% Metropolis Acceptance rule
    GAMMA_dlpi = log( mvnpdf(d, G*mu_mlpi, G*C_mlpi*G' + C_d ) ) ;
    GAMMA_dlpi_previous = log( mvnpdf(d, G*mu_mlpi_previous, G*C_mlpi_previous*G' + C_d ) ) ;
    %%% Accept or reject the  facies proposal
    if rand<exp(GAMMA_dlpi-GAMMA_dlpi_previous)
        
        facies_samples(:,it) = facies_proposal;
        
        OPERATOR = C_mlpi*G'*inv(G*C_mlpi*G' + C_d);
        mu_mlpid = mu_mlpi + OPERATOR*(d - G*mu_mlpi);        
        C_mlpid = C_mlpi - OPERATOR*G*C_mlpi;
        C_samples(:,it) = diag(C_mlpid );
        
        
        % If you want to take a sample of the continuous properties:
        model_samples(:,it) = mvnrnd(mu_mlpid, C_mlpid)';
        % else, just define it as the mean
        %model_samples(:,it) = mu_mlpid;
        
        sy_seisimic(:,it) = G*model_samples(:,it);
        
        
        mu_mlpi_previous = mu_mlpi;
        C_mlpi_previous = C_mlpi;
        log_likelihood(it) = -0.5*(d-sy_seisimic(:,it))'*inv(C_d)*(d-sy_seisimic(:,it));
        
        
    else
        facies_samples(:,it) = facies_samples(:,it-1);
        model_samples(:,it) = model_samples(:,it-1);
        sy_seisimic(:,it) = sy_seisimic(:,it-1);
        log_likelihood(it) = log_likelihood(it-1);
    end
    
end

%% Write results in the INVERSION structure
C_samples = reshape(C_samples,I,3,n_it);
INVERSION.VP.Csamples(:,:) = C_samples(:,1,it_conv:end);
INVERSION.VS.Csamples(:,:) = C_samples(:,2,it_conv:end);
INVERSION.RHOB.Csamples(:,:) = C_samples(:,3,it_conv:end);

model_samples = reshape(model_samples,I,3,n_it);
INVERSION.VP.MUsamples(:,:) = model_samples(:,1,it_conv:end);
INVERSION.VS.MUsamples(:,:) = model_samples(:,2,it_conv:end);
INVERSION.RHOB.MUsamples(:,:) = model_samples(:,3,it_conv:end);
INVERSION.VP.mean = mean(INVERSION.VP.MUsamples(:,it_conv:end)')';
INVERSION.VS.mean = mean(INVERSION.VS.MUsamples(:,it_conv:end)')';
INVERSION.RHOB.mean = mean(INVERSION.RHOB.MUsamples(:,it_conv:end)')';

[INVERSION.VP.probability,INVERSION.VP.axis,INVERSION.VP.map] = calculate_posterior_probability(INVERSION.VP.MUsamples,INVERSION.VP.Csamples,1);
[INVERSION.VS.probability,INVERSION.VS.axis,INVERSION.VS.map] = calculate_posterior_probability(INVERSION.VS.MUsamples,INVERSION.VS.Csamples,1);
[INVERSION.RHOB.probability,INVERSION.RHOB.axis,INVERSION.RHOB.map] = calculate_posterior_probability(INVERSION.RHOB.MUsamples,INVERSION.RHOB.Csamples,1);


INVERSION.FACIES.samples = facies_samples;
facies_samples_chain = facies_samples(:,it_conv:end);
for facie = 1:length(PRIOR)
    indicator = zeros(size(facies_samples_chain(:)));
    indicator(facies_samples_chain==facie) = 1;
    indicator = reshape(indicator,size(facies_samples_chain));
    facies_prob(:,facie) = sum(indicator,2)./size(indicator,2);
end
INVERSION.FACIES.prob = facies_prob;
[ ~, INVERSION.FACIES.likely] = max(facies_prob,[],2);


INVERSION.SYNTHETIC = sy_seisimic(:,2:end);

INVERSION.erro = sum( (repmat(d,1,size(INVERSION.SYNTHETIC,2)) - INVERSION.SYNTHETIC).^2,1);
INVERSION.log_likelihood = log_likelihood;


INVERSION.processing_time = toc;

end



function [mu, C] = construct_prior_facies(PRIOR, facies_sample,L)
I = length(facies_sample);
corr_12 = zeros(I,1);
corr_23 = zeros(I,1);
corr_13 = zeros(I,1);
mask_byfacies = zeros(I,I);
var_ = zeros(3*I,1);
A = covariance_matrix_exp(ones(I,1),L,1);
for facies=1:length(PRIOR)
    nu = covariance2correlation(PRIOR(facies).C);
    index = find(facies_sample==facies);
    mu(index,1) = PRIOR(facies).MU(1);
    mu(index+I,1) = PRIOR(facies).MU(2);
    mu(index+2*I,1) = PRIOR(facies).MU(3);
    var_(index) = sqrt(PRIOR(facies).C(1,1));
    var_(index+I) = sqrt(PRIOR(facies).C(2,2));
    var_(index+2*I) = sqrt(PRIOR(facies).C(3,3));
    corr_12(index) = sqrt(nu(1,2));
    corr_23(index) = sqrt(nu(2,3));
    corr_13(index) = sqrt(nu(1,3));
    mask_byfacies(index,index) = 1;
end

mask_byblock = ones(I,I);
abs_diff = abs(diff(facies_sample));
boundary_positions = find(abs_diff>0);
for block = 1:length(boundary_positions)
    mask_byblock(1:boundary_positions(block),boundary_positions(block)+1:end) = 0;
end
mask_byblock = mask_byblock.*mask_byblock';


GAMMA = [ A diag(corr_12)*A*diag(corr_12) diag(corr_13)*A*diag(corr_13);
    diag(corr_12)*A*diag(corr_12) A diag(corr_23)*A*diag(corr_23);
    diag(corr_13)*A*diag(corr_13) diag(corr_23)*A*diag(corr_23)  A ];
C = diag(var_)*GAMMA*diag(var_);

% Without the mask, the properties from different facies will be
% spaially correlated, which is not correct
% Option 1: Allow spatial correlation within the same facies
    %mask = mask_byfacies;
% Option 2: Allow spatial correlation within the same facies and within the
% same "block/layer"
    mask = mask_byblock;

C = C.*[mask  mask  mask;
    mask  mask  mask;
    mask  mask  mask];

end

function [mu, C] = construct_prior_trend_correlated_v2(PRIOR, trendVp, trendVs, trendRhob, facies_sample, L)
I = length(facies_sample);
mu = zeros(3*I,1);
var_ = zeros(3*I,1);
corr_12 = zeros(I,1);
corr_23 = zeros(I,1);
corr_13 = zeros(I,1);
mask_byfacies = zeros(I,I);
A = covariance_matrix_exp(ones(I,1),L,1);
for facies=1:length(PRIOR)
    nu = covariance2correlation(PRIOR(facies).C);
    index = find(facies_sample==facies);
    mu(index,1) = trendVp(index,facies);
    mu(index+I,1) =  trendVs(index,facies);
    mu(index+2*I,1) = trendRhob(index,facies);
    var_(index) = sqrt(PRIOR(facies).C(1,1));
    var_(index+I) = sqrt(PRIOR(facies).C(2,2));
    var_(index+2*I) = sqrt(PRIOR(facies).C(3,3));
    corr_12(index) = sqrt(abs(nu(1,2)));
    corr_23(index) = sqrt(abs(nu(2,3)));
    corr_13(index) = sqrt(abs(nu(1,3)));
    mask_byfacies(index,index) = 1;
end
GAMMA = [ A diag(corr_12)*A*diag(corr_12) diag(corr_13)*A*diag(corr_13);
    diag(corr_12)*A*diag(corr_12) A diag(corr_23)*A*diag(corr_23);
    diag(corr_13)*A*diag(corr_13) diag(corr_23)*A*diag(corr_23)  A ];
C = diag(var_)*GAMMA*diag(var_);
mask = mask_byfacies;
C = C.*[mask  mask  mask;
    mask  mask  mask;
    mask  mask  mask];
end





