function [Sigma_DKL] = calc_time_series_irr(p_j, p_vis)
%CALC_TIME_SERIES_IRR computes the time-series irreversibility of a 
%stationary masked Markov chain from jump probabilties p_j and steady-state
%occupancies p_vis of visible states 
%
% INPUTS:    
%     p_j: 2 x 2 x n_max matrix of jump probabilties
%   p_vis: matrix of visible probabilities
%
% OUTPUTS:  
%   Sigma_DKL: time-series irreversibility estimate
%
% author:   JEhrich
% version:  1.2 (2021-07-02)
% changes:  updated calculation to include visible probabilities in
% logratio

% maximum number of recorded jumps
n_max = size(p_j,3);

% reshape jump probabilities
p_12 = reshape(p_j(1,2,:),n_max,1);
p_21 = reshape(p_j(2,1,:),n_max,1);

% compute log-ratio of jump probabilities
Sigma_m = log(p_12*p_vis(2)./p_21/p_vis(1));

% delete divison-by-zero and inf entries
Sigma_m(isnan(Sigma_m) | isinf(Sigma_m)) = 0;

% sum up
Sigma_DKL = sum((p_12*p_vis(2) - p_21*p_vis(1)).*Sigma_m);

end