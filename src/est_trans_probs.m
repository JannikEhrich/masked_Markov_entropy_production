function [A, p] = est_trans_probs(x_traj)
%EST_TRANS_PROBS estimates transition probabilities and stationary
%occupation probabilities for a trajectory x_traj of a stationary Markov
%chain
%
% INPUTS:    
%   x_traj: stationary trajectory of Markov chain
%
% OUTPUTS:  
%   A: estimated transition probability matrix
%   p: stationary occupation probabilities
%
% author:   JEhrich
% version:  1.0 (2021-03-01)
% changes:  -

% size of statespace
K = length(unique(x_traj));

% stationary probabilities
p = zeros(K,1);
% count occurences
for ii = 1:K
    p(ii) = sum(x_traj == ii);
end
% normalize
p = p/sum(p);

% transition probabilities
A = zeros(K,K);
% count transitions
for ii = 1:K
    for jj = 1:K
        A(ii,jj) = sum(x_traj(2:end)==ii & x_traj(1:end-1)==jj);
    end
end
% normalize
A = A./sum(A);

end

