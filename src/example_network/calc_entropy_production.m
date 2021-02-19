function [Sigma] = calc_entropy_production(A)
%CALC_ENTROPY_PRODUCTION computes the statioanry entropy production rate of
%a the Markov chain given by transition matrix A
%
% INPUTS:    
%    A: transition probability matrix
%
% OUTPUTS:  
%   Sigma: entropy production rate
%
% author:   JEhrich
% version:  1.0 (2021-02-19)
% changes:  -

% compute steady-state probability
p = calc_steady_state(A);

% compute matrix of transition probability ratios
Sigma_m = (A.*p').*log(A./A');
% delete divison-by-zero entries
Sigma_m(isnan(Sigma_m)) = 0;
% sum up
Sigma = sum(Sigma_m(:));

end
