function [Sigma] = calc_entropy_production(A,p)
%CALC_ENTROPY_PRODUCTION computes the statioanry entropy production rate of
%a the Markov chain given by transition matrix A and steady-state
%probabilities p
%
% INPUTS:    
%    A: transition probability matrix
%    p: steady-state probability vector
%
% OUTPUTS:  
%   Sigma: entropy production rate
%
% author:   JEhrich
% version:  1.1 (2021-03-01)
% changes:  removed analytical calculation of steady-state probabilities

% compute matrix of transition probability ratios
Sigma_m = (A.*p').*log(A./A');
% delete divison-by-zero entries
Sigma_m(isnan(Sigma_m)) = 0;
% sum up
Sigma = sum(Sigma_m(:));

end
