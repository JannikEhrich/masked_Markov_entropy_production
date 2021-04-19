function [p] = calc_steady_state(A)
%CLAC_STEADY_STATE computes the steady-state distribution of a Markov chain
%for a given transition probability matrix A
%
% INPUTS:    
%    A: transition probability matrix 
%
% OUTPUTS:  
%   p: steady-state probability vector
%
% author:   JEhrich
% version:  1.0 (2021-02-19)
% changes:  -

% eigendecomposition of matrix
[V,D] = eig(A);
% steady-state has real eigenvector
V = real(V);
% find eigenvalue 1
[~, entry] = min(abs(diag(D) - 1));
% normalize eigenvector
p = V(:,entry)/sum(V(:,entry));




end

